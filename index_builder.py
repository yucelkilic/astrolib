"""
GaiaIndexBuilder — generates astrometry.net index files from Gaia DR3
for a given sky position.

Intended as a fallback when pre-built system indexes are absent or fail
to produce a WCS solution. Adapted from PhoPS/src/phops/astrometry.py.

Typical usage
-------------
from astrolib.index_builder import GaiaIndexBuilder
from pathlib import Path

builder = GaiaIndexBuilder(index_dir=Path("/data/gaia_indexes"))
ok = builder.ensure_index(ra=83.82, dec=-5.39)
if ok:
    # Pass index_dir to solve_field so solve-field picks up the new indexes
    ...
"""

from __future__ import annotations

import csv
import logging
import math
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


class GaiaIndexError(RuntimeError):
    """Raised when index generation fails in a recoverable way."""


@dataclass
class GaiaIndexBuilder:
    """
    Query Gaia DR3 around (ra, dec), apply proper-motion epoch propagation,
    and build astrometry.net index files via hpsplit + build-astrometry-index.

    Parameters
    ----------
    index_dir:
        Directory where generated index files (and the patch cache CSV)
        will be stored. Created automatically if absent.
    build_bin:
        Path or name of the ``build-astrometry-index`` binary.
    hpsplit_bin:
        Path or name of the ``hpsplit`` binary.
    catalog:
        Gaia TAP table name (default: ``gaiadr3.gaia_source``).
    radius:
        Query radius in degrees around the field centre.
    quad_scales:
        Astrometry.net quad-scale levels to generate for each healpix tile.
    cache_tolerance:
        Maximum angular distance (degrees) to consider a previous cache
        entry a valid hit (avoids redundant Gaia queries).
    apply_proper_motion:
        When True (default), propagate Gaia DR3 coordinates from epoch
        2016.0 to the observation time using proper motion and parallax.
        Requires the ``obs_time`` argument in :meth:`ensure_index`.
    """

    index_dir: Path
    build_bin: str = "build-astrometry-index"
    hpsplit_bin: str = "hpsplit"
    catalog: str = "gaiadr3.gaia_source"
    radius: float = 0.5
    quad_scales: list = field(default_factory=lambda: [0, 2, 4, 6])
    cache_tolerance: float = 0.1
    apply_proper_motion: bool = True

    def __post_init__(self) -> None:
        self.index_dir = Path(self.index_dir)
        self.index_dir.mkdir(parents=True, exist_ok=True)
        self._cache_file = self.index_dir / "generated_index.csv"

    # ── Public API ────────────────────────────────────────────────────────────

    def ensure_index(self, ra: float, dec: float, obs_time=None) -> bool:
        """
        Ensure astrometry.net index files exist for the given sky position.

        Parameters
        ----------
        ra, dec:
            Field centre in decimal degrees (ICRS).
        obs_time:
            Optional ``astropy.time.Time`` observation epoch.  Used for
            proper-motion propagation when *apply_proper_motion* is True.
            When omitted, epoch propagation is skipped.

        Returns
        -------
        bool
            True if indexes are ready (from cache or freshly generated),
            False if the process failed.
        """
        if self._cache_hit(ra, dec, obs_time):
            # Verify that the cached index files still exist on disk.
            existing = self._index_files_for(ra, dec)
            if existing:
                log.info(
                    "GaiaIndexBuilder: cache hit, %d index file(s) present — ra=%.4f dec=%.4f",
                    len(existing), ra, dec,
                )
                return True
            log.warning(
                "GaiaIndexBuilder: cache hit but no index files found on disk — regenerating ra=%.4f dec=%.4f",
                ra, dec,
            )

        log.info("GaiaIndexBuilder: building index — ra=%.4f dec=%.4f", ra, dec)
        try:
            self._check_binaries()
            patch_path = self._query_gaia(ra, dec, obs_time)
            count = self._build_indexes(patch_path, ra, dec, obs_time)
            if count == 0:
                raise GaiaIndexError("build-astrometry-index produced no output files.")
            log.info("GaiaIndexBuilder: %d index file(s) built successfully", count)
            self._save_cache(ra, dec, obs_time)
            return True
        except GaiaIndexError as exc:
            log.error("GaiaIndexBuilder: %s", exc)
            return False
        except Exception as exc:
            log.exception("GaiaIndexBuilder: unexpected error — %s", exc)
            return False

    # ── Dependency check ──────────────────────────────────────────────────────

    def _check_binaries(self) -> None:
        for name in (self.build_bin, self.hpsplit_bin):
            if not (Path(name).is_file() or shutil.which(name)):
                raise GaiaIndexError(
                    f"Required binary '{name}' not found. "
                    "Install astrometry.net and ensure it is on PATH."
                )

    # ── Cache helpers ─────────────────────────────────────────────────────────

    def _load_cache(self) -> list[tuple[float, float, float]]:
        if not self._cache_file.exists():
            return []
        entries: list[tuple[float, float, float]] = []
        with self._cache_file.open(newline="") as fh:
            for row in csv.reader(fh):
                if not row:
                    continue
                ra_c, dec_c = float(row[0]), float(row[1])
                epoch_c = float(row[2]) if len(row) >= 3 else math.nan
                entries.append((ra_c, dec_c, epoch_c))
        return entries

    def _save_cache(self, ra: float, dec: float, obs_time) -> None:
        epoch = self._epoch_bucket(obs_time) if obs_time is not None else math.nan
        with self._cache_file.open("a", newline="") as fh:
            csv.writer(fh).writerow([ra, dec, epoch])

    def _cache_hit(self, ra: float, dec: float, obs_time) -> bool:
        epoch_req = self._epoch_bucket(obs_time) if obs_time is not None else None
        for ra_c, dec_c, epoch_c in self._load_cache():
            dist = math.sqrt((ra - ra_c) ** 2 + (dec - dec_c) ** 2)
            if dist >= self.cache_tolerance:
                continue
            if epoch_req is None or math.isnan(epoch_c) or epoch_c == epoch_req:
                return True
        return False

    @staticmethod
    def _epoch_bucket(obs_time, precision: int = 1) -> float:
        return round(float(obs_time.jyear), precision)

    # ── Gaia query ────────────────────────────────────────────────────────────

    def _query_gaia(self, ra: float, dec: float, obs_time) -> Path:
        try:
            from astroquery.gaia import Gaia
        except ImportError as exc:
            raise GaiaIndexError("astroquery is required for Gaia queries") from exc

        log.info(
            "GaiaIndexBuilder: querying Gaia — ra=%.4f dec=%.4f radius=%.2f deg",
            ra, dec, self.radius,
        )
        query = f"""
            SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, parallax
            FROM {self.catalog}
            WHERE CONTAINS(
                POINT('ICRS', ra, dec),
                CIRCLE('ICRS', {ra}, {dec}, {self.radius})
            ) = 1
            ORDER BY phot_g_mean_mag
        """
        job = Gaia.launch_job_async(query, verbose=False)
        data = job.get_results()
        log.info("GaiaIndexBuilder: %d sources retrieved", len(data))

        if len(data) == 0:
            raise GaiaIndexError(
                f"Gaia query returned no stars for ra={ra:.4f} dec={dec:.4f}"
            )

        if self.apply_proper_motion and obs_time is not None:
            data = self._propagate_epoch(data, obs_time)

        suffix = f"{int(ra)}_{int(dec)}"
        patch_path = self.index_dir / f"gaia_patch_{suffix}.fits"
        data.write(str(patch_path), overwrite=True)
        log.info("GaiaIndexBuilder: patch written — %s (%d stars)", patch_path.name, len(data))
        return patch_path

    def _propagate_epoch(self, data, obs_time):
        """
        Propagate Gaia DR3 coordinates (epoch 2016.0) to the observation
        time using proper motion and parallax.  Stars with missing or
        invalid astrometric data are filtered out before propagation.
        """
        try:
            import astropy.units as u
            from astropy.coordinates import SkyCoord
            from astropy.time import Time
        except ImportError as exc:
            raise GaiaIndexError("astropy is required for epoch propagation") from exc

        def _valid_mask(column) -> np.ndarray:
            raw = getattr(column, "mask", None)
            if raw is None:
                return np.ones(len(column), dtype=bool)
            return ~np.asarray(raw, dtype=bool)

        mask = (
            _valid_mask(data["pmra"])
            & _valid_mask(data["pmdec"])
            & _valid_mask(data["parallax"])
        )
        data = data[mask]

        parallax = np.asarray(data["parallax"], dtype=float)
        valid = (parallax > 0.1) & np.isfinite(parallax)
        data = data[valid]
        parallax = parallax[valid]

        if len(data) == 0:
            raise GaiaIndexError(
                "No Gaia stars remained after proper-motion/parallax filtering."
            )

        gaia_epoch = Time(2016.0, format="jyear", scale="tdb")
        distance = (parallax * u.mas).to(u.pc, equivalencies=u.parallax())
        coords = SkyCoord(
            ra=np.asarray(data["ra"], dtype=float) * u.deg,
            dec=np.asarray(data["dec"], dtype=float) * u.deg,
            pm_ra_cosdec=np.asarray(data["pmra"], dtype=float) * u.mas / u.yr,
            pm_dec=np.asarray(data["pmdec"], dtype=float) * u.mas / u.yr,
            distance=distance,
            obstime=gaia_epoch,
        )
        propagated = coords.apply_space_motion(new_obstime=obs_time.tdb)
        data["ra"] = propagated.ra.deg
        data["dec"] = propagated.dec.deg
        data = data[np.isfinite(data["ra"]) & np.isfinite(data["dec"])]

        if len(data) == 0:
            raise GaiaIndexError("No Gaia stars remained after epoch propagation.")

        log.info("GaiaIndexBuilder: %d stars after epoch propagation", len(data))
        return data

    # ── Index generation ──────────────────────────────────────────────────────

    def _index_files_for(self, ra: float, dec: float) -> list:
        """Return existing index files matching this sky position's suffix."""
        suffix = f"{int(ra)}_{int(dec)}"
        return list(self.index_dir.glob(f"index-*_{suffix}.fits"))

    def _build_indexes(self, patch_path: Path, ra: float, dec: float, obs_time) -> int:
        """
        Build index files for all quad scales and healpix tiles.
        Returns the number of index files successfully created.
        """
        suffix = f"{int(ra)}_{int(dec)}"
        tiles = self._hpsplit(patch_path, suffix)
        built = 0

        for tile_path, tile_id in tiles:
            for scale in self.quad_scales:
                index_path = self.index_dir / f"index-550{scale:02d}-{tile_id}_{suffix}.fits"
                index_id = f"550{scale:02d}{tile_id}{abs(int(ra))}"
                cmd = [
                    self.build_bin,
                    "-i", str(tile_path),
                    "-s", "2",
                    "-P", str(scale),
                    "-E",
                    "-S", "phot_g_mean_mag",
                    "-o", str(index_path),
                    "-I", index_id,
                ]
                log.info("GaiaIndexBuilder: %s", " ".join(cmd))
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    log.warning(
                        "GaiaIndexBuilder: build-astrometry-index failed for %s:\n%s",
                        index_path.name,
                        result.stderr[:400],
                    )
                elif index_path.exists():
                    built += 1

        for tile_path, _ in tiles:
            tile_path.unlink(missing_ok=True)

        return built

    def _hpsplit(self, patch_path: Path, suffix: str) -> list[tuple[Path, str]]:
        """
        Split a FITS catalog patch into healpix tiles using hpsplit.
        Returns a list of (tile_path, tile_id) pairs.
        """
        hp_template = self.index_dir / f"gaia-hp%02i_{suffix}.fits"
        cmd = [self.hpsplit_bin, "-o", str(hp_template), "-n", "2", str(patch_path)]
        log.info("GaiaIndexBuilder: hpsplit — %s", patch_path.name)

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise GaiaIndexError(f"hpsplit failed: {result.stderr[:400]}")

        tiles = []
        for p in sorted(self.index_dir.glob(f"gaia-hp*_{suffix}.fits")):
            # Tile id is the healpix number between "gaia-hp" and the first "_"
            tile_id = p.stem.replace("gaia-hp", "").split("_")[0]
            tiles.append((p, tile_id))

        if not tiles:
            raise GaiaIndexError(f"hpsplit produced no tiles for {patch_path.name}")

        log.info("GaiaIndexBuilder: %d healpix tiles produced", len(tiles))
        return tiles

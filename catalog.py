import astropy.coordinates as coord
from astropy.table import Table
from astropy import units as u
import astropy.io.fits as fits
# from astropy.constants import c
from astropy import stats
from astropy.io import ascii
from astroquery.jplhorizons import Horizons
from .io import FileOps
from .astronomy import FitsOps
from .astronomy import TimeOps
from .astronomy import AstCalc
import numpy as np
import sep
from os import system
import math

from astroquery.xmatch import XMatch
from astroquery.vizier import Vizier

import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer

from sklearn import linear_model, datasets
from pylab import rcParams
import os

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
from .visuals import StarPlot


class Query:

    def gaia_query(self, ra_deg, dec_deg, rad_deg, max_mag=20,
                   max_coo_err=1,
                   max_sources=100):

        """
        Query Gaia DR1 @ VizieR using astroquery.vizier
        parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
        @param ra_deg: RA in degrees
        @type ra_dec: float
        @param dec_deg: DEC in degrees
        @type dec_deg: float
        @param max_mag: Limit G magnitude to be queried object(s)
        @type max_mag: float
        @max_coo_err: Max error of position
        @type max_coo_err: float
        @max_sources: Maximum number of sources
        @type max_sources: int
        @returns: astropy.table object
        """

        vquery = Vizier(columns=['Source', 'RA_ICRS',
                                 'DE_ICRS', 'e_RA_ICRS',
                                 'e_DE_ICRS', 'phot_g_mean_mag',
                                 'pmRA', 'pmDE',
                                 'e_pmRA', 'e_pmDE',
                                 'Epoch', 'Plx'],
                        column_filters={"phot_g_mean_mag":
                                            ("<{:f}".format(max_mag)),
                                        "e_RA_ICRS":
                                            ("<{:f}".format(max_coo_err)),
                                        "e_DE_ICRS":
                                            ("<{:f}".format(max_coo_err))},
                        row_limit=max_sources)

        field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                               unit=(u.deg, u.deg),
                               frame='icrs')
        return (vquery.query_region(field,
                                    width="{:f}d".format(rad_deg),
                                    catalog="I/337/gaia")[0])

    def match_catalog(self, file_name,
                      ra_keyword="ALPHA_J2000",
                      dec_keyword="DELTA_J2000",
                      catalogue='I/345/gaia2',
                      plot=False,
                      catalog_output=False, max_sources=None):
        """
        Basic circular aperture photometry of given images.
        Parameters
        ----------
        file_name : file object
             File name to be source extracted.
        ra_keyword : str
            RA keyword in catalogue.
            Default is "ALPHA_J2000"
        dec_keyword : str
            DEC keyword in catalogue.
            Default is "DEC_J2000"
        catalogue : vizer catalogue object
            Vizer catalogue name
            Default is 'I/345/gaia2'.
        plot: boolean
            Shell we plot the correlation?
        Returns
        -------
        'A dict object'
        Examples
        --------
        >>> from astrolib import catalog
        >>> from astrolib import astronomy
        >>> from astropy.table import Table

        >>> co = catalog.Query()
        >>> ac = astronomy.AstCalc()
        >>> c = ac.radec2wcs("21:05:15.250", "+07:52:6.734")

        >>> ra_list_in_deg = [c.ra.degree]
        >>> dec_list_in_deg = [c.dec.degree]

        >>> phot_object = Table([ra_list_in_deg, dec_list_in_deg], names=("ALPHA_J2000", "DELTA_J2000"))

        >>> co.match_catalog("file.fits",  catalogue="II/336/apass9",
            filter="Vmag",
            phot_object=phot_object,
            plot=True)
        """
        fo = FitsOps(file_name)
        ds = fo.source_extract()

        table = XMatch.query(cat1=ds,
                             cat2='vizier:{}'.format(catalogue),
                             max_distance=5 * u.arcsec, colRA1=ra_keyword,
                             colDec1=dec_keyword)

        if max_sources is not None:
            table = table[:max_sources]

        if catalog_output is True:
            cat_file = os.path.splitext(file_name)[0]
            ascii.write(table, "{}.csv".format(cat_file), format='csv', fast_writer=False)

        if plot is True:
            data = fo.hdu[0].data.astype(float)
            splt = StarPlot()

            splt.star_plot(data, table)

        return table

    def match_catalog_and_phot(self, file_name, ra_keyword="ALPHA_J2000",
                               dec_keyword="DELTA_J2000",
                               catalogue='I/345/gaia2',
                               filter=None,
                               phot_object=None,
                               phot_method=None,
                               plot=False,
                               catalog_output=False):
        """
        Basic circular aperture photometry of given images.
        Parameters
        ----------
        file_name : file object
             File name to be source extracted.
        ra_keyword : str
            RA keyword in catalogue.
            Default is "ALPHA_J2000"
        dec_keyword : str
            DEC keyword in catalogue.
            Default is "DEC_J2000"
        catalogue : vizer catalogue object
            Vizer catalogue name
            Default is 'I/345/gaia2'.
        filter : str
            Filter of image.
            Default is 'Vmag'.
        phot_object : astropy table
            Astropy table of objetcs that contains "ALPHA_J2000" and "DEC_J2000" columns.
            Default is None.
        phot_method : astropy table
            MAG_AUTO or None
            Default is None.
        plot: boolean
            Shell we plot the correlation?
        Returns
        -------
        'A dict object'
        Examples
        --------
        >>> from astrolib import catalog
        >>> from astrolib import astronomy
        >>> from astropy.table import Table

        >>> co = catalog.Query()
        >>> ac = astronomy.AstCalc()
        >>> c = ac.radec2wcs("21:05:15.250", "+07:52:6.734")

        >>> ra_list_in_deg = [c.ra.degree]
        >>> dec_list_in_deg = [c.dec.degree]

        >>> phot_object = Table([ra_list_in_deg, dec_list_in_deg], names=("ALPHA_J2000", "DELTA_J2000"))
        >>> phot_method = "MAG_AUTO"
        >>> co.match_catalog("file.fits",  catalogue="II/336/apass9",
            filter="Vmag",
            phot_object=phot_object,
            phot_method="MAG_AUTO",
            plot=True)
        """
        fo = FitsOps(file_name)
        object_name = fo.get_header("OBJECT")
        exptime = fo.get_header("EXPTIME")
        telescope = fo.get_header("TELESCOP")

        if filter is None:
            filter = fo.get_header("FILTER")

            if "T100" in telescope:
                if "u" in filter:
                    filter = "upmag"
                elif "g" in filter:
                    filter = "gpmag"
                elif "r" in filter:
                    filter = "rpmag"
                elif "i" in filter:
                    filter = "ipmag"
                elif "z" in filter:
                    filter = "zpmag"
                elif "U" in filter:
                    filter = "Umag"
                elif "B" in filter:
                    filter = "Bmag"
                elif "V" in filter:
                    filter = "Vmag"
                elif "R" in filter:
                    filter = "Rmag"
                elif "I" in filter:
                    filter = "Imag"
                else:
                    raise SystemExit('Filter not found!')

        ds = fo.source_extract()

        table = XMatch.query(cat1=ds,
                             cat2='vizier:{}'.format(catalogue),
                             max_distance=5 * u.arcsec, colRA1=ra_keyword,
                             colDec1=dec_keyword)
        print(table.colnames)

        if phot_method is None:
            phot_method = "MAG_AUTO"
            table['deltaMag'] = table[filter] - table["MAG_AUTO"]
        else:
            table['deltaMag'] = table[filter] - table[phot_method]

        mean, median, stddev = stats.sigma_clipped_stats(table['deltaMag'], sigma=2, maxiters=5)

        linear_zero_point = None
        linear_r2 = None
        ransac_r2 = None
        ransac_zero_point = None
        linear_calibrated_mag = []
        ransac_calibrated_mag = []

        if "FLUX" in phot_method:
            X = np.log10(np.asarray(table[phot_method]).reshape(-1, 1))
        else:
            X = np.asarray(table[phot_method]).reshape(-1, 1)
        y_ = np.asarray(table[filter]).reshape(-1, 1)
        imputer = SimpleImputer()
        y = imputer.fit_transform(y_)

        # Fit line using all data
        lr = linear_model.LinearRegression()
        lr.fit(X, y)

        # Robustly fit linear model with RANSAC algorithm
        ransac = linear_model.RANSACRegressor()
        ransac.fit(X, y)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)

        # Predict data of estimated models
        line_X = np.arange(X.min(), X.max())[:, np.newaxis]
        line_y = lr.predict(line_X)
        line_y_ransac = ransac.predict(line_X)

        # Compare estimated coefficients
        # linear_zero_point = lr.intercept_
        # linear_slope = lr.coef_
        ransac_slope = ransac.estimator_.coef_
        ransac_zero_point = ransac.estimator_.intercept_

        if phot_object is not None:
            phot_object[ra_keyword].unit = u.deg
            phot_object[dec_keyword].unit = u.deg

            for object in phot_object:
                c = coord.SkyCoord(object[ra_keyword], object[dec_keyword], frame="icrs", unit="deg")
                catalog = coord.SkyCoord(ds[ra_keyword], ds[dec_keyword], frame="icrs", unit="deg")
                max_sep = 1.0 * u.arcsec
                idx, d2d, d3d = c.match_to_catalog_3d(catalog)
                sep_constraint = d2d < max_sep
                to_be_calibrated_table = ds[idx]

                # linear_calibrated_mag.append(
                #     lr.predict(np.asarray(to_be_calibrated_table[phot_method]).reshape(-1, 1)))
                ransac_calibrated_mag.append(
                    ransac.predict(np.asarray(to_be_calibrated_table[phot_method]).reshape(-1, 1)))
        else:
            linear_calibrated_mag = [None]
            ransac_calibrated_mag = [None]

        if plot is True:
            lw = 2
            rcParams['figure.figsize'] = 8, 8
            plt.scatter(X[inlier_mask], y[inlier_mask], color='yellowgreen', marker='.',
                        label='Inliers')
            plt.scatter(X[outlier_mask], y[outlier_mask], color='gold', marker='.',
                        label='Outliers')
            # plt.plot(line_X, line_y, color='navy', linewidth=lw, label='Linear regressor')
            plt.plot(line_X, line_y_ransac, color='cornflowerblue', linewidth=lw,
                     label='RANSAC regressor')
            plt.legend(loc='lower right')
            plt.title('{} ({} - {}), {} seconds'.format(object_name, catalogue, filter, exptime))
            if "FLUX" in phot_method:
                plt.xlabel(phot_method + " (log)")
            else:
                plt.xlabel(phot_method)
            plt.ylabel(filter)
            # plt.xscale('log')
            plt.show()

        if catalog_output is True:
            cat_file = os.path.splitext(file_name)[0]
            np.savetxt("{}_ransac_inst_vs_cat.csv".format(cat_file),
                       np.concatenate((X[inlier_mask], y[inlier_mask]), axis=1), delimiter=",")
            ascii.write(table, "{}.csv".format(cat_file), format='csv', fast_writer=False)

        try:
            ransac_calibrated_mag_result = ransac_calibrated_mag[0][0][0]
        except TypeError:
            ransac_calibrated_mag_result = None

        return ({'table': table[phot_method, "MAGERR_AUTO", filter],
                 'astropy_zero_point': median,
                 'ransac_zero_point': ransac_zero_point[0],
                 'ransac_slope': ransac_slope[0][0],
                 'stddev': stddev,
                 'std_error': stddev / math.sqrt(len(X[inlier_mask])),
                 'ransac_calibrated_mag': ransac_calibrated_mag_result,
                 'ransac_calibrated_mag_err': "{}".format(stddev / math.sqrt(len(X[inlier_mask]))),
                 'ransac_equation': "{}X + {}".format(ransac_slope[0][0], ransac_zero_point[0]),
                 })

    def get_sso_ephem(self, name, epoch_start, epoch_end, epoch_step="1min", location="A84"):
        """Instantiate JPL query.

        Parameters
        ----------
        name : str, required
            Name, number, or designation of the object to be queried
        location : str or dict, optional
            Observer's location for ephemerides queries or center body
            name for orbital element or vector queries. Uses the same
            codes as JPL Horizons. If no location is provided, Earth's
            center is used for ephemerides queries and the Sun's
            center for elements and vectors queries. Arbitrary topocentic
            coordinates for ephemerides queries can be provided in the
            format of a dictionary. The
            dictionary has to be of the form {``'lon'``: longitude in
            deg (East positive, West negative), ``'lat'``: latitude in
            deg (North positive, South negative), ``'elevation'``:
            elevation in km above the reference ellipsoid, [``'body'``:
            Horizons body ID of the central body; optional; if this value
            is not provided it is assumed that this location is on Earth]}.
        epochs : scalar, list-like, or dictionary, optional
            Either a list of epochs in JD or MJD format or a dictionary
            defining a range of times and dates; the range dictionary has to
            be of the form {``'start'``:'YYYY-MM-DD [HH:MM:SS]',
            ``'stop'``:'YYYY-MM-DD [HH:MM:SS]', ``'step'``:'n[y|d|m|s]'}.
            Epoch timescales depend on the type of query performed: UTC for
            ephemerides queries, TDB for element queries, CT for vector queries.
            If no epochs are provided, the current time is used.
        id_type : str, optional
            Identifier type, options:
            ``'smallbody'``, ``'majorbody'`` (planets but also
            anything that is not a small body), ``'designation'``,
            ``'name'``, ``'asteroid_name'``, ``'comet_name'``,
            ``'id'`` (Horizons id number), or ``'smallbody'`` (find the
            closest match under any id_type), default: ``'smallbody'``


        Examples
        --------
            >>> from astroquery.jplhorizons import Horizons
            >>> eros = Horizons(id='433', location='568',
            ...              epochs={'start':'2017-01-01',
            ...                      'stop':'2017-02-01',
            ...                      'step':'1d'})
            >>> print(eros)  # doctest: +SKIP
            JPLHorizons instance "433"; location=568, epochs={'start': '2017-01-01', 'step': '1d', 'stop': '2017-02-01'}, id_type=smallbody
        """

        obj = Horizons(id=name, location=location,
                       epochs={'start': '{}'.format(epoch_start), 'stop': '{}'.format(epoch_end),
                               'step': '{}'.format(epoch_step)})
        eph = obj.ephemerides()
        return eph

    def find_skybot_objects(self,
                            odate,
                            ra,
                            dec,
                            radius=16,
                            time_travel=0,
                            observatory="A84"):

        """
        Seek and identify all the known solar system objects
        in a field of view of a given size.

        @param odate: Observation date.
        @type odate: date
        @param ra: RA of field center for search, format: degrees or hh:mm:ss
        @type ra: str
        @param dec: DEC of field center for search, format: degrees or hh:mm:ss
        @type dec: str
        @param radius: Radius.
        @type radius: float
        @param time_travel: Jump into time after given date (in hour).
        @type time_travel: float
        @param observatory: Observation code.
        @type observatory: str
        @return: str
        """

        while True:
            try:
                to = TimeOps()
                fo = FileOps()
                epoch = to.date2jd(odate) + time_travel / 24.0
                bashcmd = ("wget -q \"http://vo.imcce.fr/webservices/skybot/"
                           "skybotconesearch_query.php"
                           "?-ep={0}&-ra={1}&-dec={2}&-rm={3}&-output=object&"
                           "-loc={4}&-filter=120&-objFilter=120&-from="
                           "SkybotDoc&-mime=text\" -O skybot.cat").format(
                    epoch,
                    ra,
                    dec,
                    radius,
                    observatory)

                system(bashcmd)
                skyresult = fo.read_file_as_array("skybot.cat")
                system('rm -rf skybot.cat')
                if "No solar system object was found" not in str(skyresult):
                    tskyresult = Table(skyresult,
                                       names=('num',
                                              'name',
                                              'ra(h)',
                                              'dec(deg)',
                                              'class',
                                              'm_v',
                                              'err(arcsec)',
                                              'd(arcsec)'))
                    return (True, tskyresult)
                else:
                    return (False, str(skyresult))

            except:
                print("\nConnection Failed, Retrying..")
                continue
            break

    def known_mo_position(self, image_path=None,
                          ra=None,
                          dec=None,
                          odate=None,
                          radi=16,
                          max_mag=21):

        ac = AstCalc()
        if image_path:
            fo = FitsOps(image_path)
            if not odate:
                odate = fo.get_header('date-obs')
            else:
                odate = odate
            ra_dec = ac.center_finder(image_path, wcs_ref=True)
        elif not image_path and ra and dec and odate:
            co = coord.SkyCoord('{0} {1}'.format(ra, dec),
                                unit=(u.hourangle, u.deg),
                                frame='icrs')
            print('Target Coordinates:',
                  co.to_string(style='hmsdms', sep=':'),
                  'in {0} arcmin'.format(radi))
            ra_dec = [co.ra, co.dec]

        request0 = self.find_skybot_objects(odate,
                                            ra_dec[0].degree,
                                            ra_dec[1].degree,
                                            radius=radi)

        if request0[0]:
            asteroids = request0[1]
        elif request0[0] is False:
            print(request0[1])
            return False

        asteroids['ra_deg'] = coord.Angle(asteroids["ra(h)"], unit=u.hour)
        asteroids['dec_deg'] = coord.Angle(asteroids["dec(deg)"], unit=u.deg)

        return asteroids

    # This code adapted from vvv:
    # https://github.com/MichalZG/AsteroidsPhot/blob/master/starscoordinates.py

    def query_color(self,
                    ra,
                    dec,
                    radius=0.01,
                    min_mag=10,
                    max_mag=20,
                    max_sources=100):
        """

        Query NOMAD object
        @param ra: RA of field center for search, format: degrees or hh:mm:ss
        @type ra: str
        @param dec: DEC of field center for search, format: degrees or hh:mm:ss
        @type dec: str
        @param radius: Radius.
        @type radius: float
        @param min_mag: Minimum magnitude value of query.
        @type min_mag: float
        @param max_mag: Maximum magnitude value of query.
        @type max_mag: float
        @param max_sources: Maximum strs to be queried..
        @type max_sources: int
        @return: astropy.table
        """

        c = coord.SkyCoord(ra,
                           dec,
                           unit=(u.deg, u.deg),
                           frame='icrs')
        r = radius * u.deg

        vquery = Vizier(columns=['NOMAD1',
                                 'RAJ2000',
                                 'DEJ2000',
                                 'Bmag',
                                 'Vmag',
                                 'Rmag'],
                        column_filters={"Rmag":
                                            (">{:f}".format(min_mag)),
                                        "Rmag":
                                            ("<{:f}".format(max_mag))},
                        row_limit=max_sources)

        result = vquery.query_region(c, radius=r, catalog="NOMAD")[0]

        return (result)

    # sorts list of stars (puts the best for being comparise stars in the first place)
    def sort_stars(self, starstable, min_mag):

        starstable = starstable[starstable['Rmag'] > min_mag]

        starstable['Bmag'].fill_value = 0.0
        starstable['Vmag'].fill_value = 0.0
        starstable['Rmag'].fill_value = 0.0

        starstable = starstable.filled()
        starstable = starstable[starstable['Bmag'] != 0]
        starstable = starstable[starstable['Vmag'] != 0]
        starstable = starstable[starstable['Rmag'] != 0]

        starstable['b-v'] = starstable['Bmag'] - starstable['Vmag']
        starstable['v-r'] = starstable['Vmag'] - starstable['Rmag']
        starstable['sortby'] = abs(starstable['b-v'] -
                                   0.656) + abs(starstable['v-r'] - 0.4)

        if len(starstable) > 0:
            return (Table(np.sort(starstable, order=['sortby'])))
        else:
            print("No proper comparison star(s) found!")
            raise SystemExit

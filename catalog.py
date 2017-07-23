from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from .astronomy import FitsOps
import numpy as np
import sep
from .visuals import StarPlot


class Query:
    def gaia_query(self, ra_deg, dec_deg, rad_deg, max_mag=20,
                   max_sources=100):

        """
        Query Gaia DR1 @ VizieR using astroquery.vizier
        parameters: ra_deg, dec_deg, rad_deg: RA, Dec, field
        radius in degrees
        max_mag: upper limit G magnitude (optional)
        max_sources: maximum number of sources
        returns: astropy.table object
        """

        vquery = Vizier(columns=['Source', 'RA_ICRS',
                                 'DE_ICRS', 'e_RA_ICRS',
                                 'e_DE_ICRS', 'phot_g_mean_mag',
                                 'pmRA', 'pmDE',
                                 'e_pmRA', 'e_pmDE',
                                 'Epoch'],
                        column_filters={"phot_g_mean_mag":
                                        ("<{:f}".format(max_mag))},
                        row_limit=max_sources)
 
        field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                               unit=(u.deg, u.deg),
                               frame='icrs')
        return(vquery.query_region(field,
                                   width="{:f}d".format(rad_deg),
                                   catalog="I/337/gaia")[0])

    def match_catalog(self, file_name, radius=0.005, max_mag=20,
                      max_obj=30, plot=False):
        fo = FitsOps("108hecuba-001_R_affineremap.new")
        ds_radec, ds_xy = fo.detect_sources(skycoords=True, max_obj=30)
        qry = Query()
        gaia_list = []

        for i, do in enumerate(ds_radec):
            try:
                gaia_obj = qry.gaia_query(do.ra.degree,
                                          do.dec.degree,
                                          radius,
                                          max_mag,
                                          max_sources=1)
            except:
                pass

            gaia_list.append([gaia_obj[0][0],
                              ds_xy['x'][i],
                              ds_xy['y'][i],
                              gaia_obj[0][1],
                              gaia_obj[0][2],
                              gaia_obj[0][3],
                              gaia_obj[0][4],
                              gaia_obj[0][5],
                              gaia_obj[0][6],
                              gaia_obj[0][7],
                              gaia_obj[0][8],
                              gaia_obj[0][9],
                              gaia_obj[0][10],
                              ds_xy['flux'][i],
                              ds_xy['a'][i],
                              ds_xy['b'][i],
                              ds_xy['theta'][i],
                              ds_radec[i].ra.degree,
                              ds_radec[i].dec.degree])

            gaia_matched = np.asarray(gaia_list)

        if plot:
            data = fo.hdu[0].data.astype(float)
            bkg = sep.Background(data)
            # bkg_image = bkg.back()
            # bkg_rms = bkg.rms()
            data_sub = data - bkg
            splt = StarPlot()

            splt.star_plot(data_sub, gaia_matched)

        print(len(gaia_matched))
        return(gaia_matched)





from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import match_coordinates_sky
from astropy.table import Table
from astropy import units as u
import astropy.io.fits as fits
from astropy.constants import c
from astropy.time import Time
from astropy import stats
from .io import FileOps
from .astronomy import FitsOps
from .astronomy import TimeOps
from .astronomy import AstCalc
import numpy as np
import sep
from os import system

from astroquery.xmatch import XMatch
from astroquery.skyview import SkyView
from astroquery.vizier import Vizier
from astropy.wcs import WCS

import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from sklearn.preprocessing import Imputer


from sklearn import linear_model, datasets
from pylab import rcParams

from datetime import datetime

import os


with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   import aplpy

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
        return(vquery.query_region(field,
                                   width="{:f}d".format(rad_deg),
                                   catalog="I/337/gaia")[0])

    def match_catalog(self, file_name,
                      ra_keyword="ALPHA_J2000",
                      dec_keyword="DELTA_J2000",
                      catalogue='I/345/gaia2',
                      filter=None,
                      phot_object={None},
                      plot=False):
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

        >>> co.match_catalog("/Users/ykilic/Downloads/2059_0086_V_affineremap.fits",  catalogue="II/336/apass9",
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

        table['deltaMag'] = table[filter] - table['MAG_AUTO']

        mean, median, stddev = stats.sigma_clipped_stats(table['deltaMag'], sigma=2, maxiters=5)

        linear_zero_point = None
        linear_r2 = None
        ransac_r2 = None
        ransac_zero_point = None
        linear_calibrated_mag = []
        ransac_calibrated_mag = []

        X = np.asarray(table['MAG_AUTO']).reshape(-1, 1)
        y_ = np.asarray(table[filter]).reshape(-1, 1)
        imputer = Imputer()
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
        linear_zero_point = lr.intercept_
        linear_r2 = lr.coef_
        ransac_r2 = ransac.estimator_.coef_
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

                linear_calibrated_mag.append(
                    lr.predict(np.asarray(to_be_calibrated_table['MAG_AUTO']).reshape(-1, 1)))
                ransac_calibrated_mag.append(
                    ransac.predict(np.asarray(to_be_calibrated_table['MAG_AUTO']).reshape(-1, 1)))

        if plot is True:
            lw = 2
            rcParams['figure.figsize'] = 16, 10
            plt.scatter(X[inlier_mask], y[inlier_mask], color='yellowgreen', marker='.',
                        label='Inliers')
            plt.scatter(X[outlier_mask], y[outlier_mask], color='gold', marker='.',
                        label='Outliers')
            plt.plot(line_X, line_y, color='navy', linewidth=lw, label='Linear regressor')
            plt.plot(line_X, line_y_ransac, color='cornflowerblue', linewidth=lw,
                     label='RANSAC regressor')
            plt.legend(loc='lower right')
            plt.xlabel("MAG_AUTO")
            plt.ylabel(filter)
            plt.show()

        return ({'table': table,
                 'astropy_zero_point': median,
                 'linear_zero_point': linear_zero_point[0],
                 'ransac_zero_point': ransac_zero_point[0],
                 'linear_r2': linear_r2[0][0],
                 'ransac_r2': ransac_r2[0][0],
                 'stddev': stddev,
                 'linear_calibrated_mag': linear_calibrated_mag[0],
                 'ransac_calibrated_mag': ransac_calibrated_mag[0]})

    def match_gaia_catalog(self, file_name, radius=0.002, max_mag=20,
                           max_sources=30, plot=False):

        """
        Match detect sources with Gaia catalogue.
        @param file_name: FITS image or cat file.
        @type file_name: file or path
        @param radius: Radiys confirmation circle [in degrees]
        @type radius: float
        @param max_mag: Limit G magnitude to be queried object(s)
        @type max_mag: float
        @max_sources: Maximum number of sources
        @type max_sources: int
        @param plot: Plot detected objects?
        @type plot: boolean
        @returns: astropy.table object
        """

        fo = FitsOps(file_name)
        ds = fo.detect_sources(skycoords=True, max_sources=max_sources)

        qry = Query()
        gaia_list = []

        for i in range(len(ds)):
            try:
                gaia_obj = qry.gaia_query(ds['ra_calc'][i],
                                          ds['dec_calc'][i],
                                          radius,
                                          max_mag,
                                          max_sources=1)

                gaia_list.append([gaia_obj[0][0],
                                  ds['x'][i],
                                  ds['y'][i],
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
                                  gaia_obj[0][11],
                                  ds['mag'][i],
                                  ds['flux'][i],
                                  ds['a'][i],
                                  ds['b'][i],
                                  ds['theta'][i],
                                  ds['ra_calc'][i],
                                  ds['dec_calc'][i],
                                  (gaia_obj[0][1] -
                                   ds['ra_calc'][i]) * 3600000,
                                  (gaia_obj[0][2] -
                                   ds['dec_calc'][i]) * 3600000])

            except:
                pass

        gaia_matched = np.asarray(gaia_list)

        tgaia_matched = Table(gaia_matched, names=('id',
                                                   'x',
                                                   'y',
                                                   'ra',
                                                   'dec',
                                                   'e_ra',
                                                   'e_dec',
                                                   'g_mean_mag',
                                                   'pmra',
                                                   'pmdec',
                                                   'e_pmra',
                                                   'e_pmdec',
                                                   'epoch',
                                                   'plx',
                                                   'mag',
                                                   'flux',
                                                   'a',
                                                   'b',
                                                   'theta',
                                                   'ra_calc',
                                                   'dec_calc',
                                                   'ra_diff',
                                                   'dec_diff'))

        if plot:
            from .visuals import StarPlot
            data = fo.hdu[0].data.astype(float)
            bkg = sep.Background(data)
            data_sub = data - bkg
            splt = StarPlot()

            splt.star_plot(data_sub, tgaia_matched)

        print("Matched objects:", len(tgaia_matched))
        return (tgaia_matched)

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
                    return(True, tskyresult)
                else:
                    return(False, str(skyresult))
                
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

        return(result)

    #sorts list of stars (puts the best for being comparise stars in the first place)
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
            return(Table(np.sort(starstable, order=['sortby'])))
        else:
            print("No proper comparison star(s) found!")
            raise SystemExit
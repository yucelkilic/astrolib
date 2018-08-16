from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
from .io import FileOps
from .astronomy import FitsOps
from .astronomy import TimeOps
import numpy as np
import sep
from os import system


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

    def match_catalog(self, file_name, radius=0.002, max_mag=20,
                      max_sources=30, plot=False):

        """
        Match detect sources with Gaia catalogue.
        @param file_name: FITS image path to be search.
        @type file_name: path
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
        return(tgaia_matched)

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

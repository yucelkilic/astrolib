# -*- coding: utf-8 -*-

from astropy.io import fits
from astropy.wcs import WCS
from astropy import coordinates
from astropy import units as u
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.table import Table
from .catalog import Query
from .astronomy import FitsOps
from .astronomy import AstCalc
from .astronomy import TimeOps
import sep
import math
import numpy as np
import os


class PhotOps:

    def __init__(self, image_path):

        if image_path:
            self.image_path = image_path
            self.hdu = fits.open(self.image_path)[0]
        else:
            print("FITS image has not been provided by the user!")
            raise SystemExit

    def phot(self, x_coor, y_coor, aper_radius=3.0, gain=1.0):

        data = self.hdu.data.astype(float)

        bkg = sep.Background(data)
        # bkg_image = bkg.back()
        # bkg_rms = bkg.rms()
        data_sub = data - bkg

        flux, fluxerr, flag = sep.sum_circle(data_sub,
                                             x_coor,
                                             y_coor,
                                             aper_radius=aper_radius,
                                             err=bkg.globalrms,
                                             gain=gain)

        return({"flag": flag,
                "flux": flux,
                "fluxerr": fluxerr})

    def asteroids_phot(self, aper_radius=3.0,
                       radius=10, gain=0.57, max_mag=20):

        sb = Query()
        ac = AstCalc()
        to = TimeOps()
        fo = FitsOps(self.image_path)
        header = self.hdu.header
        w = WCS(header)

        odate = fo.get_header('date-obs')
        t1 = Time(odate.replace('T', ' '))
        exptime = fo.get_header('exptime')
        dt = TimeDelta(exptime / 2.0, format='sec')
        odate_middle = t1 + dt
        jd = to.date2jd(odate_middle.value)
        ra_dec = ac.center_finder(self.image_path, wcs_ref=True)

        request = sb.find_skybot_objects(odate_middle.value,
                                         ra_dec[0].degree,
                                         ra_dec[1].degree,
                                         radius=radius)

        if request[0]:
            asteroids = request[1]
        elif request[0] is False:
            print(request[1])
            raise SystemExit

        data = self.hdu.data.astype(float)

        bkg = sep.Background(data)
        data_sub = data - bkg

        for i in range(len(asteroids)):
            if float(asteroids['m_v'][i]) <= max_mag:
                c = coordinates.SkyCoord('{0} {1}'.format(
                    asteroids['ra(h)'][i],
                    asteroids['dec(deg)'][i]),
                                         unit=(u.hourangle, u.deg),
                                         frame='icrs')

                # asteroid's X and Y coor
                a_x, a_y = w.wcs_world2pix(c.ra.degree, c.dec.degree, 1)

                # phot asteroids
                flux, fluxerr, flag = sep.sum_circle(
                    data_sub,
                    a_x,
                    a_y,
                    aper_radius,
                    err=bkg.globalrms,
                    gain=gain)

                if flux == 0.0 or fluxerr == 0.0:
                    print("Bad asteroid selected (out of frame!)!")
                    raise SystemExit

                magt_i = ac.flux2mag(flux)
                magt_i_err = fluxerr / flux * 2.5 / math.log(10)

                min_mag_ast = float(asteroids['m_v'][i]) - 2

                if i < 1:
                    comptable = sb.query_color(c.ra.degree,
                                               c.dec.degree,
                                               5.0 / 60.0,
                                               min_mag=min_mag_ast,
                                               max_mag=19.5)

                s_comptable = sb.sort_stars(comptable, min_mag_ast)

                phot_res_list = []

                # phot comp. stars
                for j in range(len(s_comptable)):
                    # star's X and Y coor
                    s_x, s_y = w.wcs_world2pix(s_comptable['RAJ2000'][j],
                                               s_comptable['DEJ2000'][j],
                                               1)
                    # print('Circle', s_x, s_y, 10)

                    flux, fluxerr, flag = sep.sum_circle(
                        data_sub,
                        s_x,
                        s_y,
                        aper_radius,
                        err=bkg.globalrms,
                        gain=gain)

                    if flux == 0.0 or fluxerr == 0.0:
                        print("Bad star selected!")
                        raise SystemExit
                    
                    magc_i = ac.flux2mag(flux)
                    magc_i_err = fluxerr / flux * 2.5 / math.log(10)

                    magt = (float(magt_i) -
                            float(magc_i)) + s_comptable['Rmag'][j]
                    magt_err = math.sqrt(math.pow(float(magt_i_err), 2) +
                                         math.pow(float(magc_i_err), 2))

                    phot_res_list.append([asteroids['num'][i],
                                          jd,
                                          float(magt_i),
                                          float(magt_i_err),
                                          float(magc_i),
                                          float(magc_i_err),
                                          float(magt),
                                          float(magt_err),
                                          asteroids['m_v'][i],
                                          s_comptable['NOMAD1'][j],
                                          s_comptable['Rmag'][j]])

                np_phot_res = np.array(phot_res_list)

                phot_res_table = Table(np_phot_res,
                                       names=('ast_num',
                                              'jd',
                                              'magt_i',
                                              'magt_i_err',
                                              'magc_i',
                                              'magc_i_err',
                                              'magt',
                                              'magt_err',
                                              'ast_mag_cat',
                                              'nomad1',
                                              'star_Rmag'),
                                       dtype=('i4',
                                              'S25',
                                              'f8',
                                              'f8',
                                              'f8',
                                              'f8',
                                              'f8',
                                              'f8',
                                              'f8',
                                              'U20',
                                              'f8'))

                phot_res_table['magt_i'].format = '.3f'
                phot_res_table['magt_i_err'].format = '.3f'
                phot_res_table['magc_i'].format = '.3f'
                phot_res_table['magc_i_err'].format = '.3f'
                phot_res_table['magt'].format = '.3f'
                phot_res_table['magt_err'].format = '.3f'
                
                with open('{0}/{1}.txt'.format(
                        os.getcwd(),
                        asteroids['num'][i]), 'a') as f_handle:
                    f_handle.seek(0, os.SEEK_END)
                    phot_res_table.write(f_handle, format='ascii.no_header')

                print(phot_res_table)

        return(phot_res_table)

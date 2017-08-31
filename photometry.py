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
import glob
import time
import sys


class PhotOps:

    def update_progress(self, job_title, progress):
        length = 20
        block = int(round(length * progress))
        msg = "\r{0}: [{1}] {2}%".format(job_title,
                                         "#"*block + "-"*(length-block),
                                         round(progress*100, 2))
        if progress >= 1:
            msg += " DONE\r\n"

        sys.stdout.write(msg)
        sys.stdout.flush()

    def phot(self, image_path,
             x_coor, y_coor,
             aper_radius=3.0,
             gain=1.0):

        if image_path:
            hdu = fits.open(image_path)[0]
        else:
            print("FITS image has not been provided by the user!")
            raise SystemExit

        data = hdu.data.astype(float)

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

    def asteroids_phot(self, image_path, aper_radius=3.0,
                       radius=10, gain=0.57, max_mag=20):

        if ".fit" in os.path.basename(image_path):
            fitslist = sorted(glob.glob(image_path))
            if fitslist == 0:
                print('No image FITS found in the {0}'.format(image_path))
                raise SystemExit
        else:
            fitslist = sorted(glob.glob(image_path + '/*.fit?'))
            if fitslist == 0:
                print('No image FITS found in the {0}'.format(image_path))
                raise SystemExit

        for id, fitsfile in enumerate(fitslist):
            if fitsfile:
                hdu = fits.open(fitsfile)[0]
            else:
                print("FITS image has not been provided by the user!")
                raise SystemExit

            sb = Query()
            ac = AstCalc()
            to = TimeOps()
            fo = FitsOps(fitsfile)
            header = hdu.header
            w = WCS(header)

            naxis1 = fo.get_header('naxis1')
            naxis2 = fo.get_header('naxis2')
            odate = fo.get_header('date-obs')
            t1 = Time(odate.replace('T', ' '))
            exptime = fo.get_header('exptime')
            dt = TimeDelta(exptime / 2.0, format='sec')
            odate_middle = t1 + dt
            jd = to.date2jd(odate_middle.value)
            ra_dec = ac.center_finder(fitsfile, wcs_ref=True)

            request = sb.find_skybot_objects(odate_middle.value,
                                             ra_dec[0].degree,
                                             ra_dec[1].degree,
                                             radius=radius)

            if request[0]:
                asteroids = request[1]
            elif request[0] is False:
                print(request[1])
                raise SystemExit

            data = hdu.data.astype(float)

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

                    if naxis1 < a_x or naxis2 < a_y or a_x < 0 or a_y < 0:
                        continue

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

                    if i < 1 and id < 1:
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

                        if naxis1 < s_x or naxis2 < s_y or s_x < 0 or s_y < 0:
                            continue

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
                        phot_res_table.write(f_handle,
                                             format='ascii.no_header')

                    # print(phot_res_table)
            
            # Test
            time.sleep(0.1)
            self.update_progress(
                "Photometry is done for: {0}".format(fitsfile),
                id / len(fitslist))
        self.update_progress("Photometry done!", 1)
        return(True)

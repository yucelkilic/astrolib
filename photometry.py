# -*- coding: utf-8 -*-

from matplotlib import rcParams
from matplotlib.patches import Circle

from astropy.io import fits
from astropy.wcs import WCS
from astropy import coordinates
from astropy import units as u
from .catalog import Query
from .astronomy import FitsOps
from .astronomy import AstCalc
import sep
import numpy as np
import matplotlib.pyplot as plt


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
                "fluxerr":fluxerr})

    def asteroids_phot(self, aper_radius=3.0,
                       radius=10, gain=0.57, max_mag=20):

        sb = Query()
        ac = AstCalc()
        fo = FitsOps(self.image_path)
        header = self.hdu.header
        w = WCS(header)

        odate = fo.get_header('date-obs')
        ra_dec = ac.center_finder(self.image_path, wcs_ref=True)

        request = sb.find_skybot_objects(odate,
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

        phot_res_list = []

        for i in range(len(asteroids)):
            if float(asteroids['m_v'][i]) <= max_mag:
                c = coordinates.SkyCoord('{0} {1}'.format(
                    asteroids['ra(h)'][i],
                    asteroids['dec(deg)'][i]),
                                         unit=(u.hourangle, u.deg),
                                         frame='icrs')

                min_mag_ast = float(asteroids['m_v'][i]) - 2
                
                starstable = sb.query_color(c.ra.degree,
                                            c.dec.degree,
                                            10.0 / 60.0,
                                            min_mag=min_mag_ast)

                print(sb.sort_stars(starstable))

                px, py = w.wcs_world2pix(c.ra.degree, c.dec.degree, 1)
                print(asteroids['num'][i], px, py)
                
                flux, fluxerr, flag = sep.sum_circle(
                    data_sub,
                    px,
                    py,
                    aper_radius,
                    err=bkg.globalrms,
                    gain=gain)

                # magt_i = ac.flux2mag(flux)
                # magt_i_err = ac.flux2mag(fluxerr)                
        """
        sb.query_color(306.77168051, -23.47029011, 0.01, max_sources=10)
                # filename = get_pkg_data_filename(image_path)
                rcParams['figure.figsize'] = [10., 8.]
                # rcParams.update({'font.size': 10})

                data = self.hdu.data.astype(float)

                bkg = sep.Background(data)
                # bkg_image = bkg.back()
                # bkg_rms = bkg.rms()
                data_sub = data - bkg
                m, s = np.mean(data_sub), np.std(data_sub)

                ax = plt.subplot(projection=w)

                plt.imshow(data_sub, interpolation='nearest',
                           cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
                ax.coords.grid(True, color='white', ls='solid')
                ax.coords[0].set_axislabel('Galactic Longitude')
                ax.coords[1].set_axislabel('Galactic Latitude')

                overlay = ax.get_coords_overlay('icrs')
                overlay.grid(color='white', ls='dotted')
                overlay[0].set_axislabel('Right Ascension (ICRS)')
                overlay[1].set_axislabel('Declination (ICRS)')


                p = Circle((c.ra.degree, c.dec.degree), 0.001,
                           edgecolor="red",
                           facecolor='none',
                           transform=ax.get_transform('icrs'))
                ax.add_patch(p)
        # plt.gca().invert_xaxis()
        # plt.gca().invert_yaxis()
        plt.show()

        """

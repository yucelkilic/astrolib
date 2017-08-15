import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle

from .astronomy import AstCalc
from .astronomy import FitsOps
from astropy.io import fits
from astropy.table import Table
from astropy import coordinates
from astropy import units as u

import numpy as np
import sep

from astropy.wcs import WCS
# from astropy.utils.data import get_pkg_data_filename


class StarPlot:

    def star_plot(self, image_data, objects, mark_color="red"):

        """
        Source plot module.
        @param image_data: data part of the FITS image
        @type image_data: numpy array
        @param objects: Return of the detect_sources
        function with skycoords.
        @type objects: astropy.table
        @param mark_color: Color of the plot marks
        @type mark_color: str
        @returns: boolean
        """

        rcParams['figure.figsize'] = [10., 8.]
        
        # plot background-subtracted image
        fig, ax = plt.subplots()

        m, s = np.mean(image_data), np.std(image_data)
        ax.imshow(image_data, interpolation='nearest',
                  cmap='gray', vmin=m-s, vmax=m+s, origin='lower')

        # plot an ellipse for each object

        objects = Table(objects)
        
        for i in range(len(objects)):
            e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                        width=6*objects['a'][i],
                        height=6*objects['b'][i],
                        angle=objects['theta'][i] * 180. / np.pi)

            e.set_facecolor('none')
            e.set_edgecolor(mark_color)
            ax.add_artist(e)

        plt.show()

        return(True)

    def asteroids_plot(self, image_path, time_travel=1, max_mag=20.0):

        """
        Source plot module.
        @param image_data: data part of the FITS image
        @type image_data: numpy array
        @param objects: Return of the detect_sources
        function with skycoords.
        @type objects: astropy.table
        @param time_travel: Jump into time after given date (in hour).
        @type time_travel: float
        @param max_mag: Limit magnitude to be queried object(s)
        @type max_mag: float
        @param mark_color: Color of the plot marks
        @type mark_color: str
        @returns: boolean
        """

        # filename = get_pkg_data_filename(image_path)
        rcParams['figure.figsize'] = [10., 8.]
        # rcParams.update({'font.size': 10})

        hdu = fits.open(image_path)[0]
        wcs = WCS(hdu.header)

        data = hdu.data.astype(float)

        bkg = sep.Background(data)
        # bkg_image = bkg.back()
        # bkg_rms = bkg.rms()
        data_sub = data - bkg
        m, s = np.mean(data_sub), np.std(data_sub)

        ax = plt.subplot(projection=wcs)

        plt.imshow(data_sub, interpolation='nearest',
                   cmap='gray', vmin=m-s, vmax=m+s, origin='lower')
        ax.coords.grid(True, color='white', ls='solid')
        ax.coords[0].set_axislabel('Galactic Longitude')
        ax.coords[1].set_axislabel('Galactic Latitude')

        overlay = ax.get_coords_overlay('icrs')
        overlay.grid(color='white', ls='dotted')
        overlay[0].set_axislabel('Right Ascension (ICRS)')
        overlay[1].set_axislabel('Declination (ICRS)')

        fo = FitsOps(image_path)
        odate = fo.get_header('date-obs')
        
        ac = AstCalc()
        ra_dec = ac.center_finder(image_path, wcs_ref=True)

        asteroids = Table(ac.find_skybot_objects(odate,
                                                 ra_dec[0].degree,
                                                 ra_dec[1].degree,
                                                 radius=8),
                          names=('num',
                                 'name',
                                 'ra(h)',
                                 'dec(deg)',
                                 'class',
                                 'm_v',
                                 'err(arcsec)',
                                 'd(arcsec)'))

        asteroids_after = Table(ac.find_skybot_objects(
            odate,
            ra_dec[0].degree,
            ra_dec[1].degree,
            radius=8,
            time_travel=time_travel),
                                names=('num',
                                       'name',
                                       'ra(h)',
                                       'dec(deg)',
                                       'class',
                                       'm_v',
                                       'err(arcsec)',
                                       'd(arcsec)'))

        for i in range(len(asteroids)):
            if float(asteroids['m_v'][i]) <= max_mag:
                c = coordinates.SkyCoord('{0} {1}'.format(
                    asteroids['ra(h)'][i],
                    asteroids['dec(deg)'][i]),
                                         unit=(u.hourangle, u.deg),
                                         frame='icrs')

                c_after = coordinates.SkyCoord('{0} {1}'.format(
                    asteroids_after['ra(h)'][i],
                    asteroids_after['dec(deg)'][i]),
                                               unit=(u.hourangle, u.deg),
                                               frame='icrs')

                r = FancyArrowPatch(
                    (c.ra.degree, c.dec.degree),
                    (c_after.ra.degree, c_after.dec.degree),
                    arrowstyle='->',
                    mutation_scale=10,
                    transform=ax.get_transform('icrs'))

                p = Circle((c.ra.degree, c.dec.degree), 0.005,
                           edgecolor='yellow',
                           facecolor='none',
                           transform=ax.get_transform('icrs'))
                ax.text(c.ra.degree,
                        c.dec.degree - 0.007,
                        asteroids['name'][i],
                        size=12,
                        color='black',
                        ha='center',
                        va='center',
                        transform=ax.get_transform('icrs'))
                
                r.set_facecolor('none')
                r.set_edgecolor('red')
                ax.add_patch(p)
                ax.add_patch(r)
        # plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.show()
        print(asteroids)
        return(True)
        

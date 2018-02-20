import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec

from .astronomy import AstCalc
from .astronomy import FitsOps
from astropy.io import fits
from astropy.table import Table
from astropy import table
from astropy import coordinates
from astropy import units as u
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from astroquery.skyview import SkyView

import numpy as np
import sep
import os

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

    def asteroids_plot(self, image_path=None,
                       ra=None,
                       dec=None,
                       odate=None,
                       time_travel=1,
                       radius=6,
                       max_mag=20.0,
                       circle_color='yellow',
                       arrow_color='red'):

        """
        Source plot module.
        @param image_path: data part of the FITS image
        @type image_path: numpy array
        @param ra: RA coordinate of target area.
        @type ra: str in "HH MM SS"
        @param dec: DEC coordinate of target area
        @type dec: str in "+DD MM SS"
        @param radius: Radius in arcmin.
        @type radius: str
        @param odate: Ephemeris date of observation in date
        @type odate: "2017-08-15T19:50:00.95" format in str
        @param time_travel: Jump into time after given date (in hour).
        @type time_travel: float
        @param max_mag: Limit magnitude to be queried object(s)
        @type max_mag: float
        @param circle_color: Color of the asteroids marks
        @type circle_color: str
        @param arrow_color: Color of the asteroids direction marks
        @type arrow_color: str
        @returns: boolean
        """

        from .catalog import Query

        # filename = get_pkg_data_filename(image_path)
        rcParams['figure.figsize'] = [10., 8.]
        # rcParams.update({'font.size': 10})

        if image_path:
            hdu = fits.open(image_path)[0]
        elif not image_path and ra and dec and odate:
            co = coordinates.SkyCoord('{0} {1}'.format(ra, dec),
                                      unit=(u.hourangle, u.deg),
                                      frame='icrs')
            print('Target Coordinates:',
                  co.to_string(style='hmsdms', sep=':'),
                  'in {0} arcmin'.format(radius))
            try:
                server_img = SkyView.get_images(position=co,
                                                survey=['DSS'],
                                                radius=radius*u.arcmin)
                hdu = server_img[0][0]
            except:
                print("SkyView could not get the image from DSS server.")
                raise SystemExit

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

        sb = Query()
        ac = AstCalc()
        if image_path:
            fo = FitsOps(image_path)
            if not odate:
                odate = fo.get_header('date-obs')
            else:
                odate = odate
            ra_dec = ac.center_finder(image_path, wcs_ref=True)
        elif not image_path and ra and dec and odate:
            odate = odate
            ra_dec = [co.ra, co.dec]

        request0 = sb.find_skybot_objects(odate,
                                          ra_dec[0].degree,
                                          ra_dec[1].degree,
                                          radius=radius)

        if request0[0]:
            asteroids = request0[1]
        elif request0[0] is False:
            print(request0[1])
            raise SystemExit

        request1 = sb.find_skybot_objects(odate,
                                          ra_dec[0].degree,
                                          ra_dec[1].degree,
                                          radius=radius,
                                          time_travel=time_travel)

        if request1[0]:
            asteroids_after = request1[1]
        elif request1[0] is False:
            print(request1[1])
            raise SystemExit

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
                           edgecolor=circle_color,
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
                r.set_edgecolor(arrow_color)
                ax.add_patch(p)
                ax.add_patch(r)
        # plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.show()
        print(asteroids)
        return(True)

    def lc_plot_general(self, result_file_path=None,
                        xcol='jd',
                        ycol='magt_i',
                        errcol='magt_i_err',
                        mark_color="blue",
                        bar_color="red"):

        """
        Plot light curve of photometry result.
        @param result_file_path: Result file path
        @type result_file_path: path
        @param xcol: X-axis data for plotting
        @type xcol: array
        @param ycol: Y-axis data for plotting
        @type ycol: array
        @param errcol: Error bar data for plotting
        @type errcol: array
        @param mark_color: Marker color
        @type mark_color: str
        @param bar_color: Bar marker color
        @type bar_color: str
        @return: str
        """

        print("Plotting asteroid's LC...")

        fn = os.path.basename(result_file_path).split('.')[0]

        result_file = Table.read(result_file_path,
                                 format='ascii.commented_header')

        result_unique_by_jd = table.unique(result_file, keys='jd')

        rcParams['figure.figsize'] = [10., 8.]
        figlc = plt.figure(1)
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

        # Two subplots, the axes array is 1-d
        axlc1 = figlc.add_subplot(gs[0])
        axlc2 = figlc.add_subplot(gs[1])
        axlc1.set_title(fn)

        g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)

        # initialize fitters
        fit = fitting.LevMarLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,
                                                   niter=3, sigma=3.0)

        # get fitted model and filtered data
        filtered_data, or_fitted_model = or_fit(g_init,
                                                result_unique_by_jd[xcol],
                                                result_unique_by_jd[ycol])

        axlc1.plot(result_unique_by_jd[xcol],
                   filtered_data, 'ro',
                   label="filtered data")

        axlc1.errorbar(result_unique_by_jd[xcol],
                       filtered_data,
                       yerr=result_unique_by_jd[errcol],
                       fmt='o',
                       ecolor=bar_color,
                       color=mark_color,
                       capsize=5,
                       elinewidth=2)

        axlc1.invert_yaxis()
        axlc2.set_xlabel("$JD$", fontsize=12)
        axlc1.set_ylabel("$Magnitude (R)$", fontsize=12)
        axlc2.set_ylabel("$STD$", fontsize=12)

        fit = np.polyfit(result_unique_by_jd[xcol],
                         result_unique_by_jd[errcol],
                         1)
        fit_fn = np.poly1d(fit)
        axlc2.plot(result_unique_by_jd[xcol],
                   result_unique_by_jd[errcol],
                   'yo',
                   result_unique_by_jd[xcol],
                   fit_fn(result_unique_by_jd[xcol]),
                   '--k')

        axlc1.grid(True)
        axlc2.grid(True)
        axlc1.legend(loc=2, numpoints=1)

        figlc.savefig("{0}/{1}.pdf".format(os.getcwd(), fn))
        plt.show()

    def lc_plot_std_mag(self, result_file_path=None,
                        xcol='magc_i',
                        ycol='star_Rmag',
                        errcol='magc_i_err',
                        mark_color="blue",
                        bar_color="red"):

        print("Plotting asteroid's LC...")

        fn = os.path.basename(result_file_path).split('.')[0]

        result_file = Table.read(result_file_path,
                                 format='ascii.commented_header')

        result_unique_by_jd = table.unique(result_file, keys='nomad1')

        rcParams['figure.figsize'] = [10., 8.]
        figlc = plt.figure(1)
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

        # Two subplots, the axes array is 1-d
        axlc1 = figlc.add_subplot(gs[0])
        axlc2 = figlc.add_subplot(gs[1])
        axlc1.set_title(fn)

        g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
        fit = fitting.LevMarLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit,
                                                   sigma_clip,
                                                   niter=3, sigma=2.3)
        fitted_model = fit(g_init, result_unique_by_jd[xcol],
                           result_unique_by_jd[ycol])

        # get fitted model and filtered data
        filtered_data, or_fitted_model = or_fit(g_init,
                                                result_unique_by_jd[xcol],
                                                result_unique_by_jd[ycol])
        
        axlc1.plot(result_unique_by_jd[xcol], filtered_data,
                   'ro', label="filtered data")
        axlc1.plot(result_unique_by_jd[xcol], or_fitted_model,
                   'ro', label="filtered data")
        plt.show()

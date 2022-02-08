import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Ellipse
from matplotlib.patches import Circle
import matplotlib.gridspec as gridspec

from .astronomy import AstCalc
from .astronomy import FitsOps
from .io import FileOps
from astropy.io import fits
from astropy.table import Table
from astropy import table
from astropy import coordinates
from astropy import units as u
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clip, mad_std
from astroquery.skyview import SkyView
from astroquery.xmatch import XMatch
from PIL import Image

# import aplpy

import numpy as np
import sep
import os
import glob


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

        figsize = (8, 8)
        data = image_data.astype(float)
        fig, ax = plt.subplots(figsize=figsize)
        zscale = ZScaleInterval(nsamples=1000)
        ax.imshow(zscale(data), cmap="gray", aspect="auto")

        # plot an ellipse for each object
        for i in range(len(objects)):
            e = Ellipse(xy=(objects['X_IMAGE'][i], objects['Y_IMAGE'][i]),
                        width=6 * objects['A_IMAGE'][i],
                        height=6 * objects['B_IMAGE'][i],
                        angle=objects['A_IMAGE'][i] * 180. / np.pi)

            e.set_facecolor('none')
            e.set_edgecolor(mark_color)
            ax.add_artist(e)

        plt.show()

        return True

    def asteroids_plot(self,
                       image_path=None,
                       ra=None,
                       dec=None,
                       odate=None,
                       time_travel=1,
                       radi=6,
                       max_mag=20.0,
                       circle_color='yellow',
                       arrow_color='red',
                       invert_yaxis="True"):

        """
        Source plot module.
        @param image_path: data part of the FITS image
        @type image_path: numpy array
        @param ra: RA coordinate of target area.
        @type ra: str in "HH MM SS"
        @param dec: DEC coordinate of target area
        @type dec: str in "+DD MM SS"
        @param radi: Radius in arcmin.
        @type radi: float
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
        @param invert_yaxis: invert y axis or not.
        @type invert_yaxis: bool
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
                  'in {0} arcmin'.format(radi))
            try:
                server_img = SkyView.get_images(position=co,
                                                survey=['DSS'],
                                                radius=radi * u.arcmin)
                hdu = server_img[0][0]
            except Exception as e:
                print("SkyView could not get the image from DSS server.")
                print(e)
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
                   cmap='gray', vmin=m - s, vmax=m + s, origin='lower')
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
                                          radius=radi)

        if request0[0]:
            asteroids = request0[1]
        elif request0[0] is False:
            print(request0[1])
            raise SystemExit

        request1 = sb.find_skybot_objects(odate,
                                          ra_dec[0].degree,
                                          ra_dec[1].degree,
                                          radius=float(radi),
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
        if invert_yaxis == "True":
            plt.gca().invert_yaxis()
        plt.show()
        print(asteroids)
        return True

    def lc_plot_general(self,
                        result_file_path=None,
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

        result_unique_by_keys = table.unique(result_file, keys='jd')

        rcParams['figure.figsize'] = [10., 8.]
        figlc = plt.figure(1)
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

        # Two subplots, the axes array is 1-d
        axlc1 = figlc.add_subplot(gs[0])
        axlc2 = figlc.add_subplot(gs[1])
        axlc1.set_title(fn)

        filtered_data = sigma_clip(result_unique_by_keys[ycol], sigma=3,
                                   iters=10, stdfunc=mad_std)

        axlc1.errorbar(
            result_unique_by_keys[xcol][np.logical_not(filtered_data.mask)],
            result_unique_by_keys[ycol][np.logical_not(filtered_data.mask)],
            yerr=result_unique_by_keys[errcol][np.logical_not(
                filtered_data.mask)],
            fmt='o',
            ecolor=bar_color,
            color=mark_color,
            capsize=5,
            elinewidth=2)

        axlc1.invert_yaxis()
        axlc2.set_xlabel("JD", fontsize=12)
        axlc1.set_ylabel("Magnitude (R - INST)", fontsize=12)
        axlc2.set_ylabel("STD", fontsize=12)

        fit = np.polyfit(
            result_unique_by_keys[xcol][np.logical_not(filtered_data.mask)],
            result_unique_by_keys[errcol][np.logical_not(filtered_data.mask)],
            1)
        fit_fn = np.poly1d(fit)
        axlc2.plot(
            result_unique_by_keys[xcol][np.logical_not(filtered_data.mask)],
            result_unique_by_keys[errcol][np.logical_not(filtered_data.mask)],
            'yo',
            result_unique_by_keys[xcol][np.logical_not(filtered_data.mask)],
            fit_fn(result_unique_by_keys[xcol][np.logical_not(
                filtered_data.mask)]),
            '--k')

        axlc1.grid(True)
        axlc2.grid(True)
        axlc1.legend(loc=2, numpoints=1)

        figlc.savefig("{0}/{1}_jd_vs_magi_lc.pdf".format(os.getcwd(), fn))
        # plt.show()

    def lc_plot_std_mag(self, result_file_path=None,
                        xcol='magc_i',
                        ycol='star_Rmag',
                        errcol='magc_i_err',
                        mark_color="blue",
                        bar_color="red"):

        print("Plotting asteroid's LC...")

        # Fixing random state for reproducibility
        np.random.seed(19680801)

        fn = os.path.basename(result_file_path).split('.')[0]

        # Two subplots, the axes array is 1-d
        # Plotting settings
        rcParams['figure.figsize'] = [10., 8.]

        lc = plt.figure(1)
        lc_ast_std = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

        # magi vs catalogue
        lc1 = lc.add_subplot(gs[0])
        lc1.set_title(fn)
        lc1.grid(True)
        lc1.set_ylabel("Magnitude (R - NOMAD1)", fontsize=12)
        lc1.invert_yaxis()

        # magi vs STD
        lc2 = lc.add_subplot(gs[1])
        lc2.set_title(fn)
        lc2.grid(True)
        lc2.set_xlabel("Magnitude (Inst)", fontsize=12)
        lc2.set_ylabel("$STD$", fontsize=12)

        # magt vs estimated mag
        lc3 = lc_ast_std.add_subplot(gs[0])
        lc3.set_title(fn)
        lc3.legend(loc=2, numpoints=1)
        lc3.grid(True)
        lc3.invert_yaxis()
        lc3.set_xlabel("$JD$", fontsize=12)
        lc3.set_ylabel("Magnitude (R - Estimated from NOMAD1)",
                       fontsize=12)
        # Plotting settings

        result_file = Table.read(result_file_path,
                                 format='ascii.commented_header')

        # result_unique_by_keys = table.unique(result_file, keys='nomad1')
        result_unique_by_jd = table.unique(result_file, keys='jd')

        magt_std_list = []
        for jd in result_unique_by_jd['jd']:
            frame_results = result_file[(result_file['jd'] == jd)]

            # for reject outliers
            filtered_frame_results = sigma_clip(frame_results['magt_i'],
                                                sigma=3,
                                                iters=10, stdfunc=mad_std)

            # use only not rejected data (because umask used)
            filtered_f_umask = np.logical_not(filtered_frame_results.mask)

            # magci vs catalogue with error bar
            lc1.errorbar(
                frame_results[xcol][filtered_f_umask],
                frame_results[ycol][filtered_f_umask],
                yerr=frame_results[errcol][filtered_f_umask],
                fmt='o',
                ecolor=bar_color,
                color=mark_color,
                capsize=5,
                elinewidth=2)

            # magci vs catalogue fit calculation
            fit = np.polyfit(
                frame_results[xcol][filtered_f_umask],
                frame_results[ycol][filtered_f_umask],
                1)

            fit_fn = np.poly1d(fit)

            magt_to_std = fit_fn(frame_results['magt_i'][filtered_f_umask])
            magt_std_list.append([jd, magt_to_std[0], frame_results['magt_i_err'][0]])

            # magci vs catalogue fit plot
            lc1.plot(
                frame_results[xcol][filtered_f_umask],
                fit_fn(frame_results[xcol][filtered_f_umask]),
                '--k')

            # magi vs catalogue error fit calc.
            fit = np.polyfit(
                frame_results[xcol][filtered_f_umask],
                frame_results[errcol][filtered_f_umask],
                1)
            fit_fn = np.poly1d(fit)

            # magi vs STD fit plot
            lc2.plot(
                frame_results[xcol][filtered_f_umask],
                frame_results[errcol][filtered_f_umask],
                'yo',
                frame_results[xcol][filtered_f_umask],
                fit_fn(frame_results[xcol][filtered_f_umask]),
                '--k')

        # jd vs magt_std
        jd_vs_magt = np.asanyarray(magt_std_list)
        filtered_jd_vs_magt = sigma_clip(jd_vs_magt[:, 1],
                                         sigma=3,
                                         iters=10, stdfunc=mad_std)
        # use only not rejected data (because umask used)
        filtered_f_umask = np.logical_not(filtered_jd_vs_magt.mask)

        # jd vs magt plotting with error bars
        lc3.errorbar(
            jd_vs_magt[:, 0][filtered_f_umask],
            jd_vs_magt[:, 1][filtered_f_umask],
            yerr=jd_vs_magt[:, 2][filtered_f_umask],
            fmt='o',
            ecolor=bar_color,
            color=mark_color,
            capsize=5,
            elinewidth=2,
            label='{0} - R (Estimated)'.format(fn))

        lc3.legend(loc=2, numpoints=1)
        lc_ast_std.savefig("{0}/{1}_jd_vs_mag_std_lc.pdf".format(os.getcwd(), fn))

        # plt.show()

    def find_best_comp(self, result_file_path=None,
                       best_comparison_star=None):

        result_file = Table.read(result_file_path,
                                 format='ascii.commented_header')

        # read comparison star list
        # and check manual assigned comp star
        if best_comparison_star is None:
            result_unique_by_cat = table.unique(result_file, keys='nomad1')
        else:
            result_unique_by_cat = table.unique(
                result_file[(result_file['nomad1'] == best_comparison_star)],
                                                keys='nomad1')

        std_list = []
        t_c_list = []
        # calculates diff_mag for all target objects and comp. stars
        for star in result_unique_by_cat['nomad1']:
            frame_results = result_file[(result_file['nomad1'] == star)]

            # diff phot.
            frame_results['t-c'] = frame_results['magt_i'] - frame_results['magc_i']
            # error propagation
            frame_results['t-c-err'] = np.sqrt(
                np.power(frame_results['magt_i_err'], 2) + np.power(frame_results['magc_i_err'], 2))

            # extracting usefull columns
            t_c_list.append(frame_results['ast_num', 'nomad1', 'jd', 't-c', 't-c-err'])

            # calculating all t-c stars STD then adding list
            std_list.append(np.std(frame_results['t-c']))
        # calculating all STD's mean and its index number in the list
        mean_idx = (np.abs(np.asanyarray(std_list) - np.mean(std_list))).argmin()

        # choosing STD with min, mean and max stars
        diff_stats = {'min': [std_list.index(min(std_list)), min(std_list)],
                      'mean': [mean_idx, np.mean(std_list)],
                      'max': [std_list.index(max(std_list)), max(std_list)]
                      }
        # getting these diff mags and their other columns
        results = {'with_min_comp': t_c_list[diff_stats['min'][0]],
                   'with_mean_comp': t_c_list[diff_stats['mean'][0]],
                   'with_max_comp': t_c_list[diff_stats['max'][0]]
                   }

        return results

    def lc_plot_diff_mag(self, result_file_path=None,
                         best_comparison_star=None,
                         mark_color="blue",
                         bar_color="red"):

        print("Plotting asteroid's LC...")

        fn = os.path.basename(result_file_path).split('.')[0]

        # Two subplots, the axes array is 1-d
        # Plotting settings
        rcParams['figure.figsize'] = [10., 8.]

        lc_ast_diff = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

        results = self.find_best_comp(result_file_path=result_file_path,
                                      best_comparison_star=best_comparison_star)['with_mean_comp']

        filtered_jd_vs_mag_diff = sigma_clip(results['t-c'],
                                         sigma=3,
                                         iters=10, stdfunc=mad_std)

        # use only not rejected data (because umask used)
        filtered_diff_umask = np.logical_not(filtered_jd_vs_mag_diff.mask)

        # jd vs magt - magi
        lc = lc_ast_diff.add_subplot(gs[0])
        lc.set_title(fn)
        lc.legend(loc=2, numpoints=1)
        lc.grid(True)
        lc.invert_yaxis()
        lc.set_xlabel("$JD$", fontsize=12)
        lc.set_ylabel("Diff Mag. ({0} - {1})".format(
            results['ast_num'][0],
            results['nomad1'][0]),
            fontsize=12)
        # Plotting settings

        lc.errorbar(
            results['jd'][filtered_diff_umask],
            results['t-c'][filtered_diff_umask],
            yerr=results['t-c-err'][filtered_diff_umask],
            fmt='o',
            ecolor=bar_color,
            color=mark_color,
            capsize=5,
            elinewidth=2,
            label='{0} - {1}'.format(fn, results['nomad1'][0]))

        lc.legend(loc=2, numpoints=1)
        lc_ast_diff.savefig("{0}/{1}_jd_vs_diff_mag_lc.pdf".format(os.getcwd(), fn))

        # plt.show()


    def catalog_plot(self, fitsfile, catalog):

        try:
            import f2n
        except ImportError:
            print('Python cannot import f2n. Make sure f2n is installed.')
            raise SystemExit

        image = f2n.fromfits(fitsfile, verbose=False)
        image.setzscale('auto', 'auto')
        image.makepilimage('log', negative=False)

        print('\033[1;34mPlotting sources on {0}...\033[0m'.format(catalog))

        extension = os.path.splitext(os.path.basename(catalog))[1]

        if extension == '.cat':
            coordinates = np.genfromtxt(catalog, delimiter=None,
                                        comments='#')[:, [1, 2]]
        elif extension == '.txt':
            coordinates = np.genfromtxt(catalog, delimiter=None,
                                        comments='#')[:, [0, 1]]
        elif extension == '.cnd':
            coordinates = np.genfromtxt(catalog, delimiter=',', comments='#',
                                        skip_header=1)[:, [1, 2]]

        for i, coordinate in enumerate(coordinates):
            x, y = coordinate[0], coordinate[1]
            label = '{0}'.format(i + 1)
            image.drawcircle(x,
                             y,
                             r=10,
                             colour=(0, 255, 0),
                             label=label)

        image.writetitle(os.path.basename(fitsfile))
        fitshead, fitsextension = os.path.splitext(fitsfile)
        image.tonet('{0}.png'.format(fitshead))

        print('\033[1;34mAll sources plotted on: {0}.png\033[0m'.format(fitshead))

        return True

    def object_plot(self, image_path, ra, dec, mark_color="red"):

        """
        Source plot module.
        @param image_data: data part of the FITS image
        @type image_data: numpy array
        @param ra: RA coordinate of object, skycoords.
        @type ra: string
        @param dec: DEC coordinate of object, skycoords.
        @type dec: string
        @param mark_color: Color of the plot marks
        @type mark_color: str
        @returns: boolean
        """
        try:
            import f2n
        except ImportError:
            print('Python cannot import f2n. Make sure f2n is installed.')
            raise SystemExit

        if image_path:
            hdu = fits.open(image_path)[0]
        else:
            print("No image provided!")
            raise SystemExit

        wcs = WCS(hdu.header)

        # plot an ellipse for each object
        if ":" not in (ra or dec):
            co = coordinates.SkyCoord('{0} {1}'.format(ra, dec),
                                      unit=(u.deg, u.deg),
                                      frame='icrs')
        else:
            co = coordinates.SkyCoord('{0} {1}'.format(ra, dec),
                                  unit=(u.hourangle, u.deg),
                                  frame='icrs')

        print('Target Coordinates:',
              co.to_string(style='hmsdms', sep=':'))

        image = f2n.fromfits(image_path, verbose=False)
        image.setzscale('auto', 'auto')
        image.makepilimage('log', negative=False)

        ac = AstCalc()
        x, y = ac.sky2xy(image_path, ra, dec)
        label = '{0}'.format(co.to_string(style='hmsdms', sep=':'))
        image.drawcircle(x,
                         y,
                         r=10,
                         colour=(0, 255, 0),
                         label=label)

        image.writetitle(os.path.basename(image_path))
        fitshead, fitsextension = os.path.splitext(image_path)
        image.tonet('{0}.png'.format(fitshead))

        print('\033[1;34mSource plotted on: {0}.png\033[0m'.format(fitshead))

        return True

    def fits2png(self, image_path):

        """
        Source plot module.
        @param image_data: data part of the FITS image
        @type image_data: numpy array
        @param ra: RA coordinate of object, skycoords.
        @type ra: string
        @param dec: DEC coordinate of object, skycoords.
        @type dec: string
        @param mark_color: Color of the plot marks
        @type mark_color: str
        @returns: boolean
        """
        try:
            import f2n
        except ImportError:
            print('Python cannot import f2n. Make sure f2n is installed.')
            raise SystemExit

        if image_path:
            hdu = fits.open(image_path)[0]
        else:
            print("No image provided!")
            raise SystemExit

        image = f2n.fromfits(image_path, verbose=False)
        image.setzscale('auto', 'auto')
        image.makepilimage('log', negative=False)

        image.writetitle(os.path.basename(image_path))
        fitshead, fitsextension = os.path.splitext(image_path)
        image.tonet('{0}.png'.format(fitshead))

        return True

    def rota(self,
             image_path=None,
             object_name=None,
             ephemeris_file=None,
             odate=None,
             radius=None,
             srg_radius=10,
             time_travel=1,
             min_mag=0,
             max_mag=17.0,
             circle_color='yellow',
             arrow_color='red',
             invert_yaxis="True"):
        """
        Moving object trajectory plotter.
            Parameters
            ----------
            ephemeris_file: file object
                Ephemeris file.
            object_name : str
                Asteroid or moving object name.
            odate : str
                Ephemeris date of observation in date.
            min_mag : list or float
                Faintest magnitude to be plotted.
                Default is '20.0'.
            max_mag : float
                Brightest magnitude to be plotted.
                Default is '15.0'.
            circle_color : str
                Moving object mark color
            arrow_color : Trajectory color
            Returns
            -------
            'A table object or file'
        """

        from .catalog import Query

        # filename = get_pkg_data_filename(image_path)
        rcParams['figure.figsize'] = [12., 12.]
        # rcParams.update({'font.size': 10})

        fo = FileOps()
        srg = fo.srg_ephemeris_reader("/Users/ykilic/Downloads/RTT-150_20200109-1708_20200110-0423.txt")

        data_len = len(srg)
        mid = int(len(srg) / 2)
        data_mid_date = srg['Date-Time'][mid]
        ra = srg['RA2000'][mid]
        dec = srg['DECL2000'][mid]
        odate = data_mid_date

        if radius is None:
            ra_first = srg['RA2000'][0]
            dec_first = srg['DECL2000'][0]
            ra_last = srg['RA2000'][data_len-1]
            dec_last = srg['DECL2000'][data_len-1]
            c_first = coordinates.SkyCoord(ra_first, dec_first, unit=(u.hourangle, u.deg), frame='icrs')
            c_last = coordinates.SkyCoord(ra_last, dec_last, unit=(u.hourangle, u.deg), frame='icrs')
            radius = c_first.separation(c_last)
            radius = radius.arcmin

        if image_path is not None:
            hdu = fits.open(image_path)[0]
        elif (image_path is None) and ra and dec and odate:
            co = coordinates.SkyCoord('{0} {1}'.format(ra, dec),
                                      unit=(u.hourangle, u.deg),
                                      frame='icrs')
            print('Target Coordinates:',
                  co.to_string(style='hmsdms', sep=':'),
                  'in {0} arcmin'.format(radius))
            try:
                print (co)
                server_img = SkyView.get_images(position=co,
                                                survey=['DSS'],
                                                radius=radius * u.arcmin)
                hdu = server_img[0][0]
            except Exception as e:
                print("SkyView could not get the image from DSS server.")
                print(e)
                raise SystemExit

        fig = aplpy.FITSFigure(hdu, figsize=(12, 12))

        srg_c = coordinates.SkyCoord(srg['RA2000'], srg['DECL2000'], unit=(u.hourangle, u.deg), frame='icrs')
        srg['APLHA_J2000'] = srg_c.ra.degree
        srg['DELTA_J2000'] = srg_c.dec.degree

        table = XMatch.query(cat1=srg['APLHA_J2000', 'DELTA_J2000', 'RA2000', 'DECL2000'],
                             cat2='vizier:{}'.format("I/345/gaia2"),
                             max_distance= srg_radius * u.arcsec, colRA1='APLHA_J2000',
                             colDec1='DELTA_J2000')

        table_pd = table.to_pandas()

        table_pd_masked = table_pd[(table_pd['phot_g_mean_mag'] >= min_mag) &
                             (table_pd['phot_g_mean_mag'] <= max_mag)]

        # fig.show_markers(srg['APLHA_J2000'], srg['DELTA_J2000'], edgecolor='green')
        fig.show_markers(table_pd_masked['APLHA_J2000'], table_pd_masked['DELTA_J2000'], edgecolor='red')
        # fig.show_markers(srg['APLHA_J2000'], srg['DELTA_J2000'], edgecolor='red')

        srg_pd = srg.to_pandas()
        for i, element in enumerate(table_pd_masked['RA2000']):
            srg_pd = srg_pd[srg_pd.RA2000 != element]

        srg_best_positions = srg_pd

        fig.show_markers(srg_best_positions['APLHA_J2000'], srg_best_positions['DELTA_J2000'], edgecolor='blue')

        fig.show_grayscale(invert=True)
        fig.add_colorbar()
        fig.add_grid()

        fig.set_title("{} {}".format(ra, dec))

        return srg_best_positions

    def multifits2pngs(self, fitsdir):

        types = (fitsdir + '/*.fits', fitsdir + '/*.fit',
                 fitsdir + '/*.fts')  # the tuple of file types

        fits_grabbed = []
        for fits_files in types:
            fits_grabbed.extend(glob.glob(fits_files))

        if fits_grabbed:
            fits_grabbed = sorted(fits_grabbed)
        else:
            return False

        for fits_file in fits_grabbed:
            self.fits2png(fits_file)

        return True

    def make_animation(self, fitsdir):

        self.multifits2pngs(fitsdir)

        pngdir = fitsdir + '/*.png'
        png_out = fitsdir + '/animation.gif'

        img, *imgs = [Image.open(f) for f in sorted(glob.glob(pngdir))]
        img.save(fp=png_out, format='GIF', append_images=imgs,
                 save_all=True, duration=200, loop=0)

        return True

    def crop_roi(self, fits_file, source_x, source_y, roi_box=10, use_pil=False):
        body_path, ext = os.path.splitext(fits_file)

        fo = FitsOps(fits_file)

        source_roi = fo.hdu[0].data[int(source_y - roi_box):int(source_y + roi_box),
                     int(source_x - roi_box):int(source_x + roi_box)]


        norm = ImageNormalize(stretch=SqrtStretch())
        plt.axis('off')
        plt.imshow(source_roi, cmap='Greys', origin='lower', norm=norm)
        plt.savefig('{}_roi.png'.format(body_path), bbox_inches='tight', pad_inches=0, transparent=True)
        plt.close()

        return source_roi






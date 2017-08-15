import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse

from astropy.table import Table

import numpy as np


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

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Ellipse

import numpy as np


class StarPlot:

    def star_plot(self, image_data, objects, mark_color="red"):

        rcParams['figure.figsize'] = [10., 8.]
        
        # plot background-subtracted image
        fig, ax = plt.subplots()

        m, s = np.mean(image_data), np.std(image_data)
        ax.imshow(image_data, interpolation='nearest',
                  cmap='gray', vmin=m-s, vmax=m+s, origin='lower')

        # plot an ellipse for each object
        o_row, o_col = objects.shape
        
        for i in range(len(objects)):
            if o_col > 17:
                e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                            width=6*objects['a'][i],
                            height=6*objects['b'][i],
                            angle=objects['theta'][i] * 180. / np.pi)
            elif o_col == 17:
                e = Ellipse(xy=(objects[i][1], objects[i][2]),
                            width=6*objects[i][14],
                            height=6*objects[i][15],
                            angle=objects[i][16] * 180. / np.pi)

            e.set_facecolor('none')
            e.set_edgecolor(mark_color)
            ax.add_artist(e)

        plt.show()

        return(True)

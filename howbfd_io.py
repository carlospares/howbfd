from math import ceil
import matplotlib.pyplot as plt
import numpy as np

class IoManager:

    def __init__(self, plot_every, T):
        self.plot_counter = 0
        base = np.arange(0, T, plot_every)
        self.plot_times = np.append(base, T)

    def get_next_plot_time(self):
        return self.plot_times[self.plot_counter]


    def io_if_appropriate(self, x, u, t, show_plot=False, save_npy=True,
                          save_plot=True, tag=""):
        """ 
        plots (x,u) if at this timestep, t passed get_next_plot_time()
        """
        if t >= self.get_next_plot_time():
            print t
            plt.clf()
            plt.plot(x,u)
            # plt.ylim([1,4])
            plt.title(t)
            if save_plot:
                plt.savefig("{}{}.png".format(tag,t))
            if show_plot:
                plt.show()
            if save_npy:
                np.save("{}{}.npy".format(tag,t), u)
            self.plot_counter += 1

    def get_tag(self, init, perturb, equation, numflux, boundary, 
                well_balanced, N, order):
        """
        Returns a tag that summarises options used to get a simulation
        """
        return "i{}-{}e{}f{}b{}w{}n{}o{}_"\
            .format(init, perturb, equation, numflux, boundary, well_balanced,
                    N, order)

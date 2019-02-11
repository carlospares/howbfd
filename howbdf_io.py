from math import ceil
import matplotlib.pyplot as plt

def plot_if_appropriate(x, u, t, dt, plot_every=0.1, show_plots=False):
    """ 
    plots (x,u) if at this timestep, t passed an int multiple of plot_every
    plot_every=-1 never plots
    """
    if ceil((t)/plot_every) < ceil((t+dt)/plot_every):
        print t
        plt.clf()
        plt.plot(x,u)
        # plt.ylim([1,4])
        plt.title(t)
        # plt.plot(x, np.exp(x) + perturbation(x,perturb)*np.exp(t))
        plt.savefig("figs/{}.png".format(t))
        if show_plots:
            plt.show()
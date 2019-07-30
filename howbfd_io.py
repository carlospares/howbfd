# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from math import ceil
import matplotlib.pyplot as plt
import numpy as np
from initcond import InitCond
from equation import Equation
import argparse


def parse_command_line():
    """ Take all arguments from command line and return them """
    parser = argparse.ArgumentParser(description='High order well-balanced FD solver')
    parser.add_argument('-c', '--config', type=str, nargs=1, default=["howbfd_config"],
                         help='Quick configuration file')
    args = parser.parse_args()
    return args.config

def safe_name(name):
    """ Make sure the module we are going to try to import has the 
        appropriate name ("config.linear_upwind", for example) """
    if name == "howbfd_config":
        return name

    if "config" in name: # remove {... ./}config{./}
        name = name[name.rfind("config")+7:]
    if "." in name: # remove .py
        name = name[:name.rfind(".")]
    return "config." + name


class IoManager:

    def __init__(self, plot_every, T, eqn):
        self.plot_counter = 0
        base = np.arange(0, T, plot_every)
        self.plot_times = np.append(base, T)
        self.eqn = eqn

    def get_next_plot_time(self):
        return self.plot_times[self.plot_counter]


    def io_if_appropriate(self, x, u, H, t, show_plot=False, save_npy=True,
                          save_plot=True, tag=""):
        """ 
        plots (x,u) if at this timestep, t passed get_next_plot_time()
        """
        if t >= self.get_next_plot_time():
            print t
            (nvars, N) = u.shape
            plt.clf()
            self.eqn.prepare_plot(x, u, H, t)
            # for n in range(nvars):
            #     plt.plot(x,u[n,:], label="u[{}]".format(n))
            # if nvars > 1:
            # plt.legend()
            # plt.title(t)
            if save_plot:
                plt.savefig("figs/{}{}.png".format(tag,t))
            if show_plot:
                plt.show()
            if save_npy:
                np.save("npys/{}{}.npy".format(tag,t), u)
            self.plot_counter += 1

    def get_tag(self, init, perturb, equation, numflux, boundary, 
                well_balanced, N, order):
        """
        Returns a tag that summarises options used to get a simulation
        """
        return "i{}-{}e{}f{}b{}w{}n{}o{}_"\
            .format(init, perturb, equation, numflux, boundary, well_balanced,
                    N, order)

    def statistics(self, x, u, H, eqn):
        """
        Compute some statistics on u
        """
        uExact = eqn.steady(H)
        print "L1 distance to steady solution: {}".format((x[1]-x[0])*np.sum(np.abs(u-uExact),1))
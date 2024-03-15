# -*- coding: utf-8 -*-
# Carlos Parés Pulido, 2019

from math import ceil
import matplotlib.pyplot as plt
import numpy as np
from initcond import InitCond
from equation import Equation
import argparse
from parameters import Parameters

DEFAULT = "howbfd_config"
#DEFAULT = "burgers_rusanov"


def parse_command_line():
    """ Take all arguments from command line and return a Parameters object
        containing all of them """
    parser = argparse.ArgumentParser(description='High order well-balanced FD solver')
    parser.add_argument('-c', '--config', type=str, nargs=1, default=[DEFAULT],
                         help='Quick configuration file')
    parser.add_argument('-r', '--refinements', type=int, nargs=1,
                         help="If >0, duplicate mesh resolution r times and compute EOC")
    parser.add_argument('-N', type=int, nargs=1, help='Grid is of size N')
    args = parser.parse_args()
    cf = Parameters(safe_name(args.config[0]))
    
    # TO DO: repeat this for all variables controllable from command line
    if args.N is not None:
        cf.N = args.N[0]
    if args.refinements is not None:
        cf.refinements = args.refinements[0]
    
    return cf

def safe_name(name):
    """ Make sure the module we are going to try to import has the 
        appropriate name ("config.linear_upwind", for example) """
    if name == DEFAULT:
        return name
    if "config" in name: # remove {... ./}config{./}
        name = name[name.rfind("config")+7:]
    if "." in name: # remove .py
        name = name[:name.rfind(".")]
    return "config." + name


class IoManager:
    def __init__(self, eqn, cf):
        self.plot_counter = 0
        base = np.arange(0, cf.T, cf.plot_every)
        self.plot_times = np.append(base, cf.T)
        self.eqn = eqn

    def get_next_plot_time(self):
        return self.plot_times[self.plot_counter]


    def io_if_appropriate(self, x, u, H, t, cf):
        """ 
        plots (x,u) if at this timestep, t passed get_next_plot_time()
        """
        if t >= self.get_next_plot_time():
#            print t
            (nvars, N) = u.shape
            plt.clf()
            self.eqn.prepare_plot(x, u, H, t)
            dx = x[1]-x[0]
            if cf.plot_exact:
                exact = self.eqn.exact(x, t, H, cf)
                if np.any(exact): # if it's not all zeros, plot it!
#                    plt.plot(x, 1e3*(u-self.eqn.exact(x,t,cf)).T, 'r', label='error*1e3')
                    plt.plot(x, u.T,x, self.eqn.exact(x,t,H,cf).T)
                    plt.legend(['CACA'])
                    
            #plt.xlim(-1, 1)    
            #plt.ylim(-0.4, 2.5)    
            tag = self.get_tag(len(x), cf)
            if cf.save_plots:
                plt.savefig("figs/{}{}.png".format(tag,t))
            if cf.show_plots:
                plt.show()
            if cf.save_npys:
                np.save("npys/{}{}.npy".format(tag,t), u)
            self.plot_counter += 1
            
    def reset_timer(self):
        self.plot_counter = 0
                             
    def get_tag(self, N, cf):
        """
        Returns a tag that summarises options used to get a simulation
        """
        return "i{}-{}e{}H{}f{}b{}n{}o{}_"\
            .format(cf.init, cf.perturb_init, cf.equation, cf.funh, cf.nummeth, 
                    cf.boundary, N, cf.order)

    def statistics(self, x, u, H, eqn):
        """
        Compute some statistics on u.
        Probably obsolete with equation.exact() and EOC...
        """
        uSteady = eqn.steady(H)
        print "L1 distance to steady solution: {}".format((x[1]-x[0])*np.sum(np.abs(u-uSteady),1))
        
    def plot_eoc(self, dxs, errors, order=None):
        """
        Draws a loglog plot for errors wrt dxs.
        If order is not None, add a reference line of that order
        """
        plt.close() # clear possible hanging plots
        z = np.polyfit(np.log(dxs), np.log(errors), 1)
        z = z[0].round(2)
        plt.loglog(dxs, errors, 'x-', label="Best fit {}".format(z))
        if order is not None:
            plt.loglog(dxs, [ (h**order) * errors[0] / (dxs[0]**order) for h in dxs], '--', 
                        label='h^{}'.format(order))
        plt.grid()
        plt.legend()
        plt.show()

    def steady_trans_form_file(self,H,x):
        U0 = np.ones((2, len(H)))
        xx = []
        N=len(x)
        if N == 25:
            file_in = open('initial_data/analytical_sw/transcritical/trans_25.dat', 'r')
        elif N== 50:
            file_in = open('initial_data/analytical_sw/transcritical/trans_50.dat', 'r')
        elif N== 100:
            file_in = open('initial_data/analytical_sw/transcritical/trans_100.dat', 'r')
        elif N== 200:
            file_in = open('initial_data/analytical_sw/transcritical/trans_200.dat', 'r')
        elif N== 400:
            file_in = open('initial_data/analytical_sw/transcritical/trans_400.dat', 'r')
        elif N== 800:
            file_in = open('initial_data/analytical_sw/transcritical/trans_800.dat', 'r')
        elif N== 5000:
            file_in = open('initial_data/analytical_sw/transcritical/trans_5000.dat', 'r')
        for y in file_in.read().split('\n'):
            #if y.isdigit():
            xx.append(float(y))
        U0[0] = xx#2 + H
        #print U0[0]
        U0[1] = 1.53#U0[1]+24.0
        return U0


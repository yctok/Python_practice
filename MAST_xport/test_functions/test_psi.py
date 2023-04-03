# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 02:16:20 2023

@author: user
"""

# import SOLPSutils as sut
# import numpy as np
# from scipy import interpolate
# import tools as tl
# import glob

# dev = 'mast'
# shot = '027205'
# shift = 'org'
# series = ['p4_d6_9']

# basedrt, topdrt = tl.set_wdir()

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import get_cmap
from scipy import interpolate

import SOLPSutils as sut
import SOLPSxport as sxp
import inspect
try:
    import json
except:
    print('json module not available, some functionality not available')

plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.facecolor':'w'})
plt.rcParams.update({'mathtext.default': 'regular'})

# ----------------------------------------------------------------------------------------


def main(gfile_loc = None, new_filename='b2.transport.inputfile_new',
         profiles_fileloc=None, shotnum=None, ptimeid=None, prunid=None,
         nefit='tanh', tefit='tanh', ncfit='spl', chii_eq_chie = False,  # ti_eq_te = False,
         Dn_min=0.001, vrc_mag=0.0, Dn_max=200,
         chie_use_grad = False, chii_use_grad = False, new_b2xportparams = True,
         chie_min = 0.01, chii_min = 0.01, chie_max = 400, chii_max = 400,
         reduce_Ti_fileloc = None, update_old_last10s = False,
         fractional_change = 1, exp_prof_rad_shift = 0, ti_fileloc = None,
         impurity_list = [], use_existing_last10=False, plot_xport_coeffs=True,
         plotall=False, verbose=False, figblock=False,
         ti_decay_len=0.015, te_decay_len = None, ne_decay_len = None,
         ti_decay_min=1, te_decay_min = 1, ne_decay_min = 1e18, shift= 1):
    """
    Driver for the code, returns an object of class 'SOLPSxport'

    Inputs:
      rundir            Location of SOLPS run directory
                        --> depricated, assume you are in the directory (this is required for b2plot calls)
      gfile_loc         Location of gfile used for SOLPS grid (for mapping SOLPS grid to psin)
      profiles_fileloc  (optional) Location of a .pkl file with saved Tom's tools profile
                        fit coefficients (produced by getProfDBPedFit from SOLPSutils)
                        OR if not a .pkl file, can use a pfile as long as profiles extend far enough into SOL
      shotnum,ptimeid,prunid  Profile identifying shot number, timeid and runid (from Tom's tools)
                        (this is uneccessary if you have a .pkl file given in profiles_fileloc)
      xxfit             Fit function used in each profile fit (xx => ne, te and nc)
      chii_eq_chie      Set to True to ignore Ti profile and just set chi_i = chi_e (not implemented yet!!)
      Dn_min            Set a floor for the allowable particle diffusion coefficient
      Dn_max            Set a maximum for the allowable particle diffusion coefficient
      chie_max          Set a maximum for the allowable electron energy diffusion coefficient
      vrc_mag           Hard-coded carbon impurity pinch, for trying to match nC profiles
                        (leave zero unless you also change the function within calcXportCoeffs)
      ti_decay_len      Decay length (at the outboard midplane) for imposed exponential falloff
                        for experimental Ti, beginning at separatrix (impurity CER is incorrect in SOL)
      te_decay_len      ""
      ne_decay_len      ""
      ti_decay_min      far-SOL Ti to decay to (eV)
      te_decay_min      far-SOL Te to decay to (eV)
      ne_decay_min      far-SOL ne to decay to (m^-3)
      chie/i_use_grad   Use ratio of the gradients for new values of chi_e/i, rather than fluxes
      new_b2xportparams Produces updated b2.transport.parameters so that D, X are set in PFR to match first
                        radial cell of SOL (default is on)
      use_existing_last10  Set to True if you have already run 2d_profiles to produce .last10 files
                           in the run folder to save time. Otherwise this will call 2d_profiles so
                           that you don't accidentally use .last10 files from a previous iteration
                           with out-of-date SOLPS profiles
      reduce_Ti_fileloc Set to None to use T_D = T_C from MDS+ profile fit
                        *On GA clusters (Iris and Omega), example file is located here:
                        '/fusion/projects/results/solps-iter-results/wilcoxr/T_D_C_ratio.txt'
      update_old_last10s  Set to True to copy the last10 files to last10.old for comparison with the next iteration
      fractional_change Set to number smaller than 1 if the incremental change is too large and
                        you want to take a smaller step
      exp_prof_rad_shift: Apply a radial shift to experimental profiles
                        (in units of psi_n, positive shifts profiles outward so separatrix is hotter)
      ti_fileloc        Separate file with Ti data (overwrites previous fits)
      impurity_list     List of all the impurities included in the plasma simulation
                        (not tested yet for anything other than 'c')
      plot_xport_coeffs Plot the SOLPS and experimental profiles, along with the previous
                        and next iteration of transport coefficients
      plotall           Produce a bunch more plots from the subroutines used in the process
      verbose           Set to True to print a lot more outputs to terminal
      figblock          Set to True if calling from command line and you want to see figures
      shift             Set the amount of shifting the major radius
      
    Returns:
      Object of class 'SOLPSxport', which can then be used to plot, recall, or modify the saved data
      and rewrite a new b2.transport.inputfile
    """
    if 'json' in sys.modules:
        # Write dict of last call arguments as json file
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)
        argvals = {}
        for arg in args:
            argvals[arg] = values[arg]
        json.dump(argvals,open("SOLPSxport_args.last",'w'))

    if shotnum is None:
        try:
            shotnum = int(gfile_loc[gfile_loc.rfind('g')+1:gfile_loc.rfind('.')])
        except ValueError:
            if verbose: print("Can't determine shot number, setting to 0")
            shotnum = 0
    if ptimeid is None:
        try:
            ptimeid = int(gfile_loc[gfile_loc.rfind('.')+1:])
        except ValueError:
            if verbose: print("Can't determine time slice, setting to 0")
            ptimeid = 0
    ptimeid = int(ptimeid)
    shotnum = int(shotnum)

    print("Initializing SOLPSxport")
    xp = sxp.SOLPSxport(workdir=os.getcwd(), gfile_loc=gfile_loc, impurity_list=impurity_list)

    print("Reading SOLPS output")
    try:
        dsa = sut.read_dsa("dsa")
        b2mn = sut.scrape_b2mn("b2mn.dat")        
        geo = sut.read_b2fgmtry("../baserun/b2fgmtry")
        state = sut.read_b2fstate("b2fstate")
        xport = sut.read_transport_files(".", dsa=dsa, geo=geo, state=state)
    except:
        print('Failed to read output directly, will try using b2plot')
        sut.set_b2plot_dev(verbose=verbose)
        xp.b2plot_ready = True
        dsa = None
        b2mn = None
        geo = None
        state = None
        xport = None
        
    print("Running calcPsiVals")
    try:
        xp.calcPsiVals(plotit=plotall,geo=geo,b2mn=b2mn,dsa=dsa,shift=shift)
    except Exception as err:
        print('Exiting from SOLPSxport_dr\n')
        sys.exit(err)
    
    return xp
    
if __name__ == '__main__':
    import argparse

    py3_9 = (sys.version_info[0] >= 3 and sys.version_info[1] >= 9)

    parser = argparse.ArgumentParser(description='Generate new b2.transport.inputfile files for SOLPS',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-g', '--gfileloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
    parser.add_argument('-p', '--profilesloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
    parser.add_argument('-s', '--shotnum', help='shot number; default = None', type=str, default=None)
    parser.add_argument('-t', '--timeid', help='time of profile run; default = None', type=str, default=None)
    parser.add_argument('-r', '--runid', help='profile run id; default = None', type=str, default=None)
    parser.add_argument('-i', '--tiratiofile', help='File location for Ti/TD ratio; default = None', type=str, default=None)
    parser.add_argument('-d', '--tdfileloc', help='File location for TD; default = None', type=str, default=None)
    parser.add_argument('-f', '--fractional_change', help='Fractional change to transport coefficients; default = 1',
                        type=float, default=1)
    parser.add_argument('-sh', '--shift', help='shift of major radius; default = 1', type=float, default=1)
    if py3_9:
        parser.add_argument('--chii_eq_chie', action='store_true', default=False)
        # parser.set_defaults(chii_eq_chie=False)
        parser.add_argument('--chie_use_grad', action='store_true', default=False)
        # parser.set_defaults(chie_use_grad=False)
        parser.add_argument('--chii_use_grad', action='store_true', default=False)
        # parser.set_defaults(chii_use_grad=True)

    args = parser.parse_args()

    if not py3_9:
        args.chii_eq_chie = False
        args.chie_use_grad = False
        args.chii_use_grad = False

    _ = main(gfile_loc=args.gfileloc, profiles_fileloc=args.profilesloc,
             shotnum=args.shotnum, ptimeid=args.timeid, prunid=args.runid,
             ti_fileloc=args.tdfileloc,
             chii_eq_chie=args.chii_eq_chie, chie_use_grad=args.chie_use_grad, chii_use_grad=args.chii_use_grad,
             reduce_Ti_fileloc=args.tiratiofile, fractional_change=args.fractional_change, figblock=True, shift=args.shift)
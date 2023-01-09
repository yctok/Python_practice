# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 18:50:02 2023

@author: user
"""

"""
This driver is used to generate new b2.transport.inputfile files for SOLPS
that will come closer to matching the experimental upstream profiles.
If things go well, after several iterations, this will produce radial
ne, Te and Ti profiles in SOLPS at the outboard midplane that match the
experimental profiles provided.
It does this through the generation of an object of class 'SOLPSxport'

The routine "main" runs the sequence in the correct order and returns
the SOLPSxport object for additional plotting and debugging if necessary


Instructions for command line call:

->Source SOLPS-ITER setup file for b2plot calls
->Navigate to run directory
$ python ~/Pytools/SOLPSxport_dr.py -g <gfile location> -s <shot number> -t <profile time ID> -r <profile run ID>
or if you already have a saved profiles file (.pkl):
$ python ~/Pytools/SOLPSxport_dr.py -g <gfile location> -p <profiles file location>


Requirements:
-Code depends on the SOLPSxport class contained in SOLPSxport.py, which in turn
depends on routines in SOLPSutils.py
-An existing SOLPS run, which will be used to estimate the next iteration's
radial transport coefficients
-Experimental Te, ne and Ti profiles, preferably saved in the MDSplus database ("Tom's tools")
-Source setup.ksh from a SOLPS-ITER distribution that can run b2plot before launching
ipython and importing this module (b2plot is used it grab SOLPS data)

For best results, use core density boundary condition matched to experimental data
at innermost SOLPS grid location

This usually requires several iterations to match the experimental profiles well

Once you've figured out all of the settings you need and are iterating on a run to
converge to the solution for transport coefficients, the routine "increment_run" can
be useful to do this all quickly

R.S. Wilcox, J.M. Canik and J.D. Lore 2020
contact: wilcoxrs@ornl.gov

Reference for this procedure:
https://doi.org/10.1016/j.jnucmat.2010.11.084
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import get_cmap

import SOLPSxport as sxp

plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.facecolor':'w'})
plt.rcParams.update({'mathtext.default': 'regular'})



import argparse

parser = argparse.ArgumentParser(description='Generate new b2.transport.inputfile files for SOLPS',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-g', '--gfileloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
parser.add_argument('-p', '--profilesloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
parser.add_argument('-s', '--shotnum', help='shot number; default = None', type=str, default=None)
parser.add_argument('-t', '--timeid', help='time of profile run; default = None', type=str, default=None)
parser.add_argument('-r', '--runid', help='profile run id; default = None', type=str, default=None)
parser.add_argument('-i', '--tifileloc', help='File location for Ti/TD ratio; default = None', type=str, default=None)
parser.add_argument('-d', '--device', help='Device being simulated with SOLPS; default = None', type=str, default=None)
args = parser.parse_args()


   
def main_MAMmod(gfile_loc = None, profiles_fileloc=None, shotnum=None, 
                ptimeid=None, prunid=None, device='MAST'):

 """
 Driver for the code using Osborne profile fits saved in MDSplus

 Inputs:
   rundir            Location of SOLPS run directory
                     --> depricated, assume you are in the directory (this is required for b2plot calls)
   gfile_loc         Location of gfile used for SOLPS grid (for mapping SOLPS grid to psin)
   profiles_fileloc  (optional) Location of a .pkl file with saved Tom's tools profile
                     fit coefficients (produced by getProfDBPedFit from SOLPSutils)
   shotnum,ptimeid,prunid  Profile identifying shot number, timeid and runid (from Tom's tools)
                     (this is uneccessary if you have a .pkl file given in profiles_fileloc)
   xxfit             Fit function used in each profile fit (xx => ne, te and nc)
   Dn_min            Set a floor for the allowable particle diffusion coefficient
   Dn_max            Set a maximum for the allowable particle diffusion coefficient
   chie_max          Set a maximum for the allowable electron energy diffusion coefficient
   vrc_mag           Hard-coded carbon impurity pinch, for trying to match nC profiles
                     (leave zero unless you also change the function within calcXportCoeffs)
   ti_decay_len      Decay length (at the outboard midplane) for imposed exponential falloff
                     for experimental Ti, beginning at separatrix
                     (since we know Ti measurement from CER is incorrect in SOL)
   ke/i_use_grad     Use ratio of the gradients for new values of chi_e/i, rather than fluxes
   use_existing_last10  Set to True if you have already run 2d_profiles to produce .last10 files
                        in the run folder to save time. Otherwise this will call 2d_profiles so
                        that you don't accidentally use .last10 files from a previous iteration
                        with out-of-date SOLPS profiles
   reduce_Ti_fileloc Set to None to use T_D = T_C from MDS+ profile fit
   carbon            Set to False if SOLPS run includes D only
                     note: this routine is not yet generalized to anything other than D or D+C
   plot_xport_coeffs Plot the SOLPS and experimental profiles, along with the previous
                     and next iteration of transport coefficients
   plotall           Produce a bunch more plots from the subroutines used in the process
   verbose           Set to True to print a lot more outputs to terminal
   figblock          Set to True if calling from command time and you want to see figures

 Returns:
   Object of class 'SOLPSxport', which can then be called to plot or recall the saved data
 """  
    new_filename='b2.transport.inputfile_new'
    nefit='tanh'
    tefit='tanh'
    ncfit='spl'
    Dn_min=0.001
    vrc_mag=0.0
    ti_decay_len=0.015
    Dn_max=20,
    ke_use_grad = False
    ki_use_grad = False
    chie_min = 0.01
    chii_min = 0.01
    chie_max = 200
    chii_max = 200,
    reduce_Ti_fileloc='/fusion/projects/results/solps-iter-results/wilcoxr/T_D_C_ratio.txt',
    carbon=True
    use_existing_last10=False
    plot_xport_coeffs=True
    plotall=False
    verbose=False
    figblock=False
    device='diii-d'
    computer='iris'
    TieqTe=False
    TT=True





    if shotnum is None: 
        shotnum = int(gfile_loc[-12:-6])
    else: 
        shotnum = int(shotnum)

    if device == 'MAST':
        carbon = False
        TieqTe = True
        reduce_Ti_fileloc = None
        ti_decay_len = None
        TT = False


    print("Initializing SOLPSxport")
    xp = sxp.SOLPSxport(workdir=os.getcwd(), gfile_loc=gfile_loc, carbon_bool=carbon)
    print("Running calcPsiVals")
    xp.calcPsiVals(plotit=plotall)
    print("Running getSOLPSlast10Profs")
    xp.getSOLPSlast10Profs(plotit=plotall, use_existing_last10=use_existing_last10)
    xp.loadProfDBPedFit_MAMmod(profiles_fileloc, shotnum, ptimeid, prunid, TT=TT, verbose=True)
    print("Populating PedFits")
    # xp.loadPopulateCMODPedFits(profiles_fileloc, shotnum, npsi=250, plotit=plotall)
    xp.populatePedFits_MAMmod(nemod=nefit, temod=tefit, ncmod=ncfit, npsi=250, TT=TT, TieqTe=TieqTe, plotit=plotall)
    print("Getting flux profiles")

    if carbon:
        print("Running getSOLPSCarbonProfs")
        xp.getSOLPSCarbonProfs(plotit=plotall)
    xp.getSOLPSfluxProfs(plotit=plotall)

    print("Running calcXportCoeff")
    xp.calcXportCoef(plotit=plotall or plot_xport_coeffs, reduce_Ti_fileloc=reduce_Ti_fileloc, Dn_min=Dn_min,
                     ti_decay_len=ti_decay_len, TieqTe=TieqTe, TT=TT, vrc_mag=vrc_mag, verbose=verbose, Dn_max=Dn_max,
                     chii_min=chii_min, chii_max=chii_max, chie_min=chie_min, chie_max=chie_max, figblock=figblock)
    
    print("Running writeXport")
    xp.writeXport(new_filename=new_filename, ke_use_grad=ke_use_grad, ki_use_grad=ki_use_grad)

    return xp


if __name__ == '__main__':
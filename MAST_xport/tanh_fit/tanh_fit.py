# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:27:35 2023

@author: user
"""

import numpy as np
import argparse
import mtanh_fitting as mf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate

parser = argparse.ArgumentParser(description='load MAST data file')
parser.add_argument('-d','--mastfile_loc', type=str, default=None, help='location MAST data file')
args = parser.parse_args()

def read_mastfile(mastfile_loc):
    with open(mastfile_loc, mode='r') as dfile:
        lines = dfile.readlines()
    
    profiles = {}
    nlines_tot = len(lines)
    psi_n = np.zeros(nlines_tot)
    ne = np.zeros(nlines_tot)
    te = np.zeros(nlines_tot)
    i = 0
    
    while i < nlines_tot:
        r_line = lines[i].split()
        psi_n[i] = float(r_line[0])
        ne[i] = float(r_line[1])*pow(10, -20)
        te[i] = float(r_line[3])/1000
        i += 1

    profiles['psi_normal'] = psi_n
    profiles['electron_density(10^20/m^3)'] = ne
    profiles['electron_temperature(KeV)'] = te
    return profiles

if __name__ == '__main__':
    mast_dat_dict = read_mastfile(args.mastfile_loc)
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    te = mast_dat_dict['electron_temperature(KeV)']
    
    ne_fit = mf.super_fit_osbourne(psi, ne, maxfev=400)
    ne_superfit = mf.best_osbourne(psi, ne, plot=True)
    
    plt.figure(1)
    plt.plot(psi, ne_fit[0], color='r', label= 'electron density fit')
    plt.scatter(psi, ne, label= 'electron density experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron density: ${n_e}$ (m$^{-3}$)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron density',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
    
    plt.show()
    
    # plt.figure(2)
    # plt.plot(psi, te_fit[0], color='r', label= 'electron density fit')
    # plt.scatter(psi, te, label= 'electron density experiment data')
    
    # plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    # plt.ylabel('Electron temperature: ${T_e}$ (eV)', fontdict={"family":"Times New Roman","size": 20})
    # plt.title('Electron temperature',fontdict={"family":"Times New Roman","size": 20})
    # plt.legend()

    
    
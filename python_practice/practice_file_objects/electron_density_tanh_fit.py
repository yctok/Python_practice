# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:33:01 2023

@author: Yi-Cheng
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

# def TANH(r, r0, h, d, b, m):
#     return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

#[0,3.5e20,0.005,1e18,1e21]

# def tanh_f(r, b, h, r0, d):
#     return b+ h*(np.tanh((r0-r)/d))

def tanh_close(x, h, a, m, b, a_2, m_2, b_2):
    return h+ (a_1)*(np.tanh((m_1)*x+ b_1)) +  m*(r0-r-d)*np.heaviside(r0-r-d, 1)

# x_data = np.linspace(0, 1.5, 100)
# tanh_fit = tanh_f(x_data, 0.2, 0.25, 0.976, 0.1)

# plt.plot(x_data, tanh_fit, color='r')
# plt.show()

p0 = [0.2, 0.2, -33.3, 33.2, 0.05, -33.3, 200]

if __name__ == '__main__':
    mast_dat_dict = read_mastfile(args.mastfile_loc)
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    # x_data = np.linspace(0, 1.2, 100)
    # tanh_fit = tanh_close(x_data, 0.2, 0.2, -33.3, 33.2, 0.05, -33.3, 200)        
    popt, pcov = curve_fit(tanh_close, psi, ne, p0)
    print(popt)    
    tanh_cfit = tanh_close(psi, popt[0], popt[1], popt[2], popt[3],popt[4], popt[5], popt[6])
    # plt.plot(x_data, tanh_fit, color='r')
    plt.plot(psi, tanh_cfit, color='r')
    plt.scatter(psi, ne)
    plt.show()
    
    

# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 19:47:42 2023

@author: user
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


def tanh_fit2(x, h, a_1, m_1, b_1, a_2, m_2, b_2):
    return h+ (a_1)*(np.tanh((m_1)*x+ b_1)) + (a_2)*(np.tanh((m_2)*x+ b_2))


p0 = [0.35, 0.05, -45.3, 45.2, 0.5, -2.2, 1.32]

if __name__ == '__main__':
    mast_dat_dict = read_mastfile(args.mastfile_loc)
    psi = mast_dat_dict['psi_normal']
    te = mast_dat_dict['electron_temperature(KeV)']
    # x_data = np.linspace(0.2, 1.1, 500)
    # temp_fit = tanh_fit2(x_data, 0.35, 0.05, -45.3, 45.2, 0.5, -2.2, 1.32)        
    popt, pcov = curve_fit(tanh_fit2, psi, te, p0)
    print(popt)  
    x_model = np.linspace(min(psi), max(psi), 300)
    tanh_cfit = tanh_fit2(x_model, popt[0], popt[1], popt[2], popt[3],popt[4], popt[5], popt[6])
    # plt.plot(x_data, temp_fit, color='r')
    plt.plot(x_model, tanh_cfit, color='r')
    plt.scatter(psi, te)
    plt.show()
    
    

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 00:02:41 2023

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

def tanh(r, r0, h, d, b, m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)


n_tot = 200
p0 = [0.97, 0.6, 0.01, -0.3, 3/14]
p1 = [0.95, 0.2, 0.02, 0, 6/7]

if __name__ == '__main__':
    mast_dat_dict = read_mastfile(args.mastfile_loc)
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    te = mast_dat_dict['electron_temperature(KeV)']
    popt_ne, pcov_ne = curve_fit(tanh, psi, ne, p0)
    print(popt_ne)
    popt_te, pcov_te = curve_fit(tanh, psi, te, p1)
    print(popt_te) 
    x_model = np.linspace(min(psi), max(psi), n_tot)
    tanh_ne_fit = tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
    tanh_te_fit = tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
    
    plt.figure(1)
    plt.plot(x_model, tanh_ne_fit, color='r')
    plt.scatter(psi, ne)
    
    plt.figure(2)
    plt.plot(x_model, tanh_te_fit, color='r')
    plt.scatter(psi, te)
       
    plt.show()
    
    w_datalist = []   
    for j in range(n_tot):
        w_list =[]
        w_list.append("{: .6f}".format(x_model[j]))
        w_list.append("{: .6f}".format(tanh_ne_fit[j]))
        w_list.append("{: .6f}".format(tanh_te_fit[j]))
        w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
        w_datalist.append(w_writelist)
   
    with open('fit_27205_275.dat', 'w') as f:
        for l,w_line in enumerate(w_datalist):   
            f.writelines(w_line + "\n")
    
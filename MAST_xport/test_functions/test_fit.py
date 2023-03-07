# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 14:16:24 2023

@author: user
"""

import numpy as np
import argparse
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

def tanh_fit(x, h, a, m, b):
    return h+ a*(np.tanh(m*x+ b))

def tanh_fit2(x, h, a_1, m_1, b_1, a_2, m_2, b_2):
    return h+ (a_1)*(np.tanh((m_1)*x+ b_1)) + (a_2)*(np.tanh((m_2)*x+ b_2))

def tanh_fit3(x, h, a_1, m_1, b_1, a_2, m_2, b_2, a_3, m_3, b_3):
    return h+ (a_1)*(np.tanh((m_1)*x+ b_1)) + (a_2)*(np.tanh((m_2)*x+ b_2)) + (a_3)*(np.tanh((m_3)*x+ b_3))

n_tot = 1000
p0 = [0.25, 0.131, -93.3, 92.67, 0.0385, -15.2, 10.25]
p1 = [0.35, 0.05, -45.3, 45.2, 0.5, -2.2, 1.32]
psi_solps_high = 1.09145489
psi_solps_low = 0.56591402

if __name__ == '__main__':
    mast_dat_dict = read_mastfile(args.mastfile_loc)
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    te = mast_dat_dict['electron_temperature(KeV)']

    popt_ne, pcov_ne = curve_fit(tanh_fit2, psi, ne, p0)
    with np.printoptions(precision=3, suppress=True):
        print(popt_ne)

    popt_te, pcov_te = curve_fit(tanh_fit2, psi, te, p1)
    with np.printoptions(precision=3, suppress=True):
        print(popt_te)

    x_model = np.linspace(min(psi), psi_solps_high, n_tot)
    # x_model = np.linspace(min(psi), psi_solps_high, n_tot)
    tanh_ne_fit = tanh_fit2(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4], popt_ne[5], popt_ne[6])
    tanh_te_fit = tanh_fit2(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4], popt_te[5], popt_te[6])
    
    ne_f = interpolate.interp1d(psi, ne, kind= 'linear', fill_value = 'extrapolate')
    y_ne = ne_f(x_model)
    
    error_ne = np.zeros(n_tot)
    for k in range(n_tot):
        error_ne[k] = abs(tanh_ne_fit[k] - y_ne[k])*100/y_ne[k]
        
    te_f = interpolate.interp1d(psi, te, kind= 'linear', fill_value = 'extrapolate')
    y_te = te_f(x_model)
    
    error_te = np.zeros(n_tot)
    for m in range(n_tot):
        error_te[m] = abs(tanh_te_fit[m] - y_te[m])*100/y_te[m]
    
    
    plt.figure(1)
    plt.plot(x_model, error_ne, color='r')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Error persentage: (%)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron density fit error',fontdict={"family":"Times New Roman","size": 20})
    
    
    plt.figure(2)
    plt.plot(x_model, error_te, color='r')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Error persentage: (%)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron temperature fit error',fontdict={"family":"Times New Roman","size": 20})
       
    
    plt.figure(3)
    plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
    plt.scatter(psi, ne, label= 'electron density experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron density: ${n_e}$ (m$^{-3}$)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron density',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
    
    
    plt.figure(4)
    plt.plot(x_model, tanh_te_fit, color='r', label= 'electron temperature fit')
    plt.scatter(psi, te, label= 'electron temperature experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron temperature: ${T_e}$ (eV)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron temperature',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
     
    plt.show()
    
    fitfiles = {}
    psi_n = np.zeros(n_tot)
    ne = np.zeros(n_tot)
    te = np.zeros(n_tot)
    i = 0
    
    while i < n_tot:      
        psi_n[i] = x_model[i]
        ne[i] = tanh_ne_fit[i]*pow(10, 20)
        te[i] = tanh_te_fit[i]*1000
        i += 1

    fitfiles['psi_normal'] = psi_n
    fitfiles['electron_density(m^3)'] = ne
    fitfiles['electron_temperature(eV)'] = te
    
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

        

    
    
    
    

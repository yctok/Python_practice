# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 00:02:41 2023

@author: user
"""

import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tools as tl
import glob

# parser = argparse.ArgumentParser(description='load MAST data file')
# parser.add_argument('-d','--mastfile_loc', type=str, default=None, help='location MAST data file')
# args = parser.parse_args()

# def read_mastfile(mastfile_loc):
#     with open(mastfile_loc, mode='r') as dfile:
#         lines = dfile.readlines()
    
#     profiles = {}
#     nlines_tot = len(lines)
#     psi_n = np.zeros(nlines_tot)
#     ne = np.zeros(nlines_tot)
#     te = np.zeros(nlines_tot)
#     i = 0
    
#     while i < nlines_tot:
#         r_line = lines[i].split()
#         psi_n[i] = float(r_line[0])
#         ne[i] = float(r_line[1])*pow(10, -20)
#         te[i] = float(r_line[3])/1000
#         i += 1

#     profiles['psi_normal'] = psi_n
#     profiles['electron_density(10^20/m^3)'] = ne
#     profiles['electron_temperature(KeV)'] = te
#     return profiles

# def tanh(r, r0, h, d, b, m):
#     return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)


basedrt, topdrt = tl.set_wdir()
dev = 'mast'
shot = '027205'

g_loc = glob.glob('{}/{}/{}/yag_*'.format(topdrt, dev, shot))
print(g_loc)

n_tot = 200
p0 = [0.97, 0.6, 0.01, -0.3, 3/14]
p1 = [0.95, 0.2, 0.02, 0, 6/7]

if __name__ == '__main__':
    mast_dat_dict = tl.read_mastfile(g_loc[-1])
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    te = mast_dat_dict['electron_temperature(KeV)']
    popt_ne, pcov_ne = curve_fit(tl.tanh, psi, ne, p0)
    print(popt_ne)
    popt_te, pcov_te = curve_fit(tl.tanh, psi, te, p1)
    print(popt_te) 
    x_model = np.linspace(min(psi), max(psi), n_tot)
    tanh_ne_fit = tl.tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
    tanh_te_fit = tl.tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
    
    gnexp = np.gradient(tanh_ne_fit)
    
    
    plt.figure(1)
    plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
    plt.scatter(psi, ne, label= 'electron density experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron density',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
    
    
    
    plt.figure(2)
    plt.plot(x_model, tanh_te_fit, color='r', label= 'electron temperature fit')
    plt.scatter(psi, te, label= 'electron temperature experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron temperature: ${T_e}$ (KeV)', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron temperature',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
    
    
    plt.figure(3)
    plt.plot(x_model, gnexp, color='r', label= 'gradiant')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Times New Roman","size": 20})
    plt.ylabel('Electron density gradiant', fontdict={"family":"Times New Roman","size": 20})
    plt.title('Electron density gradiant',fontdict={"family":"Times New Roman","size": 20})
    plt.legend()
    
    plt.show()
    
    w_datalist = []
    filename = 'fit_027205_275.dat'
    fdir = '{}/{}/{}'.format(basedrt, dev, filename)
    for j in range(n_tot):
        w_list =[]
        w_list.append("{: .6f}".format(x_model[j]))
        w_list.append("{: .6f}".format(tanh_ne_fit[j]))
        w_list.append("{: .6f}".format(tanh_te_fit[j]))
        w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
        w_datalist.append(w_writelist)
   
    with open(fdir, 'w') as f:
        for l,w_line in enumerate(w_datalist):   
            f.writelines(w_line + "\n")
    
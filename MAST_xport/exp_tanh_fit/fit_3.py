# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 00:02:41 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tools as tl
import glob

basedrt, topdrt = tl.set_wdir()
dev = 'mast'
shot = '027205'

g_loc = glob.glob('{}/{}/{}/yag_*'.format(topdrt, dev, shot))
print(g_loc)

n_tot = 200
# solps_psi = 1.09145489
solps_psi = 1.1
p0 = [0.97, 0.6, 0.01, -5, 3/14]
p1 = [0.95, 0.2, 0.02, -5, 6/7]

if __name__ == '__main__':
    mast_dat_dict = tl.read_mastfile(g_loc[-1])
    psi = mast_dat_dict['psi_normal']
    ne = mast_dat_dict['electron_density(10^20/m^3)']
    ne_er = mast_dat_dict['density error(10^20/m^3)']
    te = mast_dat_dict['electron_temperature(KeV)']
    te_er = mast_dat_dict['temperature error(10^20/m^3)']
    
    popt_ne, pcov_ne = curve_fit(tl.tanh, psi, ne, p0)
    print(popt_ne)
    popt_te, pcov_te = curve_fit(tl.tanh, psi, te, p1)
    print(popt_te)
      
    x_model = np.linspace(min(psi), max(psi), n_tot)
    tanh_ne_fit = tl.tanh(x_model, popt_ne[0], popt_ne[1], popt_ne[2], popt_ne[3], popt_ne[4])
    tanh_te_fit = tl.tanh(x_model, popt_te[0], popt_te[1], popt_te[2], popt_te[3], popt_te[4])
    
    shift = 0.006
    psi_sh = psi + shift
    
    sh_opt_ne, sh_cov_ne = curve_fit(tl.tanh, psi_sh, ne, p0)
    print(sh_opt_ne)
    sh_opt_te, sh_cov_te = curve_fit(tl.tanh, psi_sh, te, p1)
    print(sh_opt_te)
    
    
    x_sh = np.linspace(min(psi_sh), max(psi_sh), n_tot)
    sh_ne_fit = tl.tanh(x_sh, sh_opt_ne[0], sh_opt_ne[1], sh_opt_ne[2], sh_opt_ne[3], sh_opt_ne[4])
    sh_te_fit = tl.tanh(x_sh, sh_opt_te[0], sh_opt_te[1], sh_opt_te[2], sh_opt_te[3], sh_opt_te[4])
    
    
    gnexp = np.gradient(tanh_ne_fit)
    
    plt.figure(1)
    plt.plot(x_model, tanh_ne_fit, color='r', label= 'electron density fit')
    plt.errorbar(psi, ne, ne_er,fmt="o", label= 'electron density experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Calibri","size": 20})
    plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)', fontdict={"family":"Calibri","size": 20})
    plt.title('Electron density',fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    
    
    plt.figure(2)
    plt.plot(x_model, tanh_te_fit, color='r', label= 'electron temperature fit')
    plt.errorbar(psi, te, te_er, fmt= "o", label= 'electron temperature experiment data')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Calibri","size": 20})
    plt.ylabel('Electron temperature: ${T_e}$ (KeV)', fontdict={"family":"Calibri","size": 20})
    plt.title('Electron temperature',fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    plt.figure(3)
    plt.plot(x_sh, sh_ne_fit,'-o', color='r', label= 'electron density fit with shift')
    plt.plot(x_model, tanh_ne_fit,'-o', color='b', label= 'electron density fit')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Calibri","size": 20})
    plt.ylabel('Electron density: ${n_e}$ (10$^{20}$*m$^{-3}$)', fontdict={"family":"Calibri","size": 20})
    plt.title('Electron density',fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    plt.figure(4)
    plt.plot(x_sh, sh_te_fit,'-o', color='r', label= 'electron temperature fit with shift')
    plt.plot(x_model, tanh_te_fit,'-o', color='b', label= 'electron temperature fit')
    
    plt.xlabel('Magnetic flux coordinate: ${\psi_N}$', fontdict={"family":"Calibri","size": 20})
    plt.ylabel('Electron temperature: ${T_e}$ (KeV)', fontdict={"family":"Calibri","size": 20})
    plt.title('Electron temperature',fontdict={"family":"Calibri","size": 20})
    plt.legend()
    
    
    plt.show()
    
    # w_datalist = []
    # filename = 'fit_027205_275.dat'
    # fdir = '{}/{}/{}'.format(basedrt, dev, filename)
    # for j in range(n_tot):
    #     w_list =[]
    #     w_list.append("{: .6f}".format(x_model[j]))
    #     w_list.append("{: .6f}".format(tanh_ne_fit[j]))
    #     w_list.append("{: .6f}".format(tanh_te_fit[j]))
    #     w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
    #     w_datalist.append(w_writelist)
   
    # with open(fdir, 'w') as f:
    #     for l,w_line in enumerate(w_datalist):   
    #         f.writelines(w_line + "\n")
    
    
    w_datalist = []
    filename = 'wsh_027205_275.dat'
    fdir = '{}/{}/{}'.format(basedrt, dev, filename)
    for j in range(n_tot):
        w_list =[]
        w_list.append("{: .6f}".format(x_sh[j]))
        w_list.append("{: .6f}".format(sh_ne_fit[j]))
        w_list.append("{: .6f}".format(sh_te_fit[j]))
        w_writelist = ' '.join(str(y)+ "\t" for y in w_list)
        w_datalist.append(w_writelist)
   
    with open(fdir, 'w') as f:
        for l,w_line in enumerate(w_datalist):   
            f.writelines(w_line + "\n")
    
    
    
    
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 23:33:53 2023

@author: user
"""

import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

psi_solps = [0.56591402, 0.59635553, 0.65365526, 0.70507622, 0.75083496,
       0.79132874, 0.82698182, 0.85806206, 0.88490369, 0.90789432,
       0.92738632, 0.94367313, 0.95706941, 0.96795829, 0.97677538,
       0.9838775 , 0.98955578, 0.99415907, 0.99803803, 1.002408  ,
       1.00753157, 1.01263476, 1.01772166, 1.02279374, 1.02785249,
       1.03288158, 1.03794617, 1.04306613, 1.04817989, 1.05328886,
       1.05838546, 1.06347049, 1.06855367, 1.07363646, 1.07872032,
       1.08380671, 1.08889011, 1.09145489]

dsa = [-0.10817856, -0.09990042, -0.08461333, -0.07136643, -0.05987293,
       -0.04989893, -0.04124898, -0.03377955, -0.02737025, -0.02191068,
       -0.01729219, -0.01343244, -0.01025504, -0.00766848, -0.00557757,
       -0.00389362, -0.00254028, -0.00144284, -0.00052025,  0.00052025,
        0.00174129,  0.00295871,  0.00417348,  0.00538548,  0.0065947 ,
        0.00779984,  0.00900968,  0.01022689,  0.01144393,  0.01266685,
        0.01389157,  0.01511209,  0.01633297,  0.01755502,  0.01877859,
        0.02000397,  0.02123034,  0.02184967]

plt.figure(1)
plt.plot(psi_solps, dsa, color='r')
plt.show()

def read_fitfile(fitfile_loc):
    with open(fitfile_loc, mode='r') as dfile:
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
        ne[i] = float(r_line[1])*pow(10, 20)
        te[i] = float(r_line[2])*1000
        i += 1

    profiles['psi_normal'] = psi_n
    profiles['electron_density(m^3)'] = ne
    profiles['electron_temperature(eV)'] = te
    return profiles


fit_dat_dict = read_fitfile('fit_27205_275.dat')
teexppsi = fit_dat_dict['psi_normal']
neexppsi = fit_dat_dict['psi_normal']
neexp = fit_dat_dict['electron_density(m^3)']
teexp = fit_dat_dict['electron_temperature(eV)']


#check electron density

psi_to_dsa_func = interpolate.interp1d(psi_solps, dsa, fill_value = 'extrapolate')
dsa_neprofile = psi_to_dsa_func(neexppsi)

plt.figure(2)
plt.plot(neexppsi, dsa_neprofile, color='r')
plt.show()

gnexp = np.gradient(neexp) / np.gradient(dsa_neprofile)

gnexp_dsafunc = interpolate.interp1d(dsa_neprofile, gnexp, kind='linear', fill_value = 'extrapolate')
gnexp_solpslocs = gnexp_dsafunc(dsa)

plt.figure(3)
plt.plot(dsa, gnexp_solpslocs, color='r')
plt.show()

expden_dsa_func = interpolate.interp1d(dsa_neprofile, neexp, kind='linear', fill_value = 'extrapolate')

ne_decay_len_end = (expden_dsa_func(dsa[-2]) - expden_dsa_func(dsa[-1])) / \
    np.mean([expden_dsa_func(dsa[-1]), expden_dsa_func(dsa[-2])])

print(ne_decay_len_end)


# Check electron temperature

psi_to_dsa_func = interpolate.interp1d(psi_solps, dsa, fill_value = 'extrapolate')
dsa_teprofile = psi_to_dsa_func(teexppsi)

plt.figure(4)
plt.plot(teexppsi, dsa_teprofile, color='r')
plt.show()

gteexp = np.gradient(teexp) / np.gradient(dsa_teprofile)

gteexp_dsafunc = interpolate.interp1d(dsa_teprofile, gteexp, kind='linear', fill_value = 'extrapolate')
gteexp_solpslocs = gteexp_dsafunc(dsa)

plt.figure(5)
plt.plot(dsa, gteexp_solpslocs, color='r')
plt.show()

expTe_dsa_func = interpolate.interp1d(dsa_teprofile, teexp, kind='linear', fill_value = 'extrapolate')

te_decay_len_end = (expTe_dsa_func(dsa[-2]) - expTe_dsa_func(dsa[-1])) / \
    np.mean([expTe_dsa_func(dsa[-1]), expTe_dsa_func(dsa[-2])])

print(te_decay_len_end)
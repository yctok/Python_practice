# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:09:12 2023

@author: Yi-Cheng
"""

"Give a try"

import numpy as np


#import pickle

# exp_file = np.genfromtxt('yag_27205_275.dat',delimiter=' ')

# psi_n = []
# exp_ne = []
# err_ne =[]
# exp_te = []
# err_te = []

# j=0

# while j< 60:
#     a = exp_file[j,0]
#     psi_n.append(a)
#     b = exp_file[j,2]
#     exp_ne.append(b)
#     c = exp_file[j,4]
#     err_ne.append(c)
#     d = exp_file[j,10]
#     exp_te.append(d)
#     e = exp_file[j,16]
#     err_te.append(e)  
#     j= j+ 1

# psi_dat = {'psi_data': psi_n}
# ne_dat = {'electron_density': exp_ne, 'electron_density_error': err_ne}
# te_dat = {'electron_temperature': exp_te, 'electron_temperature_error': err_te}

# mast_dat = {'psi_data': psi_dat, 'density_data': ne_dat, 'temperature_data': te_dat}

# with open("my_pickle.pkl", "wb") as f:
#     pickle.dump(mast_dat, f)

# with open("192012_3250_e8099.pkl","rb") as pkl_example:
#     myexample = pickle.load(pkl_example)



with open(yag_27205_275.dat, mode='r') as pfile:
    lines = pfile.readlines()

lin1 = lines[linestart].split()


# def read_pfile(pfile_loc):
#     """
#     Read in the kinetic profiles from a p file to be used as inputs (successfully tested 2018/1/3)

#     Returns a dictionary with a non-intuitive set of keys (units are included)
    
#     ** Note: pfiles don't normally go into the SOL **
#     """
#     with open(pfile_loc, mode='r') as pfile:
#         lines = pfile.readlines()

#     profiles = {}
#     nprofs = 0  # counter for total number of profiles so far
#     linestart = 0  # counter for which line to start at for each profile
#     nlines_tot = len(lines)
    
#     p_1 = 'psi_normal'
#     p_2 = 'electron_density'
#     p_3 = 'electron_temperature'

#     while True:
#         # Read the header line for each profile first
#         lin1 = lines[linestart].split()
#         npts_prof = int(lin1[0])

        

#         # Generate and populate the profile arrays
#         x = np.zeros(npts_prof)
#         y = np.zeros(npts_prof)
#         dy = np.zeros(npts_prof)
#         for i in range(npts_prof):
#             split_line = lines[linestart + i + 1].split()
#             x[i] = float(split_line[0])
#             y[i] = float(split_line[1])
#             dy[i] = float(split_line[2][:-1])


#         nprofs += 1
#         linestart += 1 + npts_prof

#         if linestart >= nlines_tot:
#             break
#     return profiles

    # Check if all psinorms are the same, consolidate if so (they are, don't bother separating)

    # condense = True
    # psinorm = None
    # for k in profiles.keys():
    #     if k is None or k=='':
    #         continue
    #
    #     if k[:4] == 'psin':
    #         if psinorm is None:
    #             psinorm = profiles[k]
    #
    #         if max(abs(profiles[k] - psinorm)) > 1e-5:
    #             condense = False
    #             break

    # if condense:
    #     profiles = {key: value for key, value in profiles.items()
    #                 if key[:4] != 'psin' or key is None or key==''}
    #     profiles['psinorm'] = psinorm






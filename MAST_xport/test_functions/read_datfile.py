# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:09:12 2023

@author: Yi-Cheng
"""

"Give a try"



# import numpy as np

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
#     profiles['electron_density'] = ne
#     profiles['electron_temperature'] = te
#     return profiles

# data = read_mastfile('yag_27205_275.dat')
# print(data['psi_normal'])

"try rfind"

"""
rfind(): return the position of the last substring found in the string
"""
# a = 'mississippi'
# a.rfind('is')
# print(a.rfind('is'))

def find_time(mastfile_loc):
    timeid = mastfile_loc[mastfile_loc.rfind('_')+1:mastfile_loc.rfind('.')]
    return timeid

print(find_time('yag_27205_275.dat'))




    
    



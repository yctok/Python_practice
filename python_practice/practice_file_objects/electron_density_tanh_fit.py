# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:33:01 2023

@author: Yi-Cheng
"""

import numpy as np
import argparse

parser = argparse.ArgumentParser(description='load mast data')
parser.add_argument('-d', '--mastfile_loc', type=str, default=None, help='location of profs_*.pkl saved profile file')
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
    result = read_mastfile(args.mastfile_loc)

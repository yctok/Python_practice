# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:59:41 2019

@author: rmreksoatmodjo

Collection of general Tools to perform oft-repeated SOLPS data analyis and post-processing tasks
"""
import os
import numpy as np

def set_wdir(): #Function to set correct Working Directory Path depending on which machine is in use
    if os.environ['OS'] == 'Windows_NT':
        if os.environ['USERNAME'] == 'Yi-Cheng':
            basedrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Simulation_Data"
            topdrt = r"C:/Users/Yi-Cheng/Documents/SOLPS_Data/Experimental_Data"
        elif os.environ['USERNAME'] == 'user':
            basedrt = r"C:/Users/user/Documents/SOLPS data/simulation data"
            topdrt = r"C:/Users/user/Documents/SOLPS data/experiment data"
    else:
        print('please add new directory in tools')
    
    return basedrt, topdrt

def a_number(text):
    name = text.split('\\',1)[1]
    nu = int(name.split('_')[2])

    return [nu, text]

def ex1_number(text):
    name = text.split('\\',1)[1]
    nu = int(name.split('_')[3])

    return [nu, text]

def new_number(text):
    name = text.split("/",-1)[-2]
    nu = int(name.split('_')[2])

    return [nu, text]

def unit_dic():
    unit = {
        'ne3da.last10': ['Electron density', 'Electron density: ${n_e}$ (m$^{-3}$)'],
        'te3da.last10': ['Electron temperature', 'Electron temperature: ${T_e}$ (eV)'],
        'an3da.last10':['Neutral density', 'Neutral density: ${n_D}$ (m$^{-3}$)']
        
        }
    return unit

def read_mastfile(mastfile_loc):
    with open(mastfile_loc, mode='r') as dfile:
        lines = dfile.readlines()
    
    profiles = {}
    nlines_tot = len(lines)
    psi_n = np.zeros(nlines_tot)
    ne = np.zeros(nlines_tot)
    ne_er = np.zeros(nlines_tot)
    te = np.zeros(nlines_tot)
    te_er = np.zeros(nlines_tot)
    
    i = 0
    
    while i < nlines_tot:
        r_line = lines[i].split()
        psi_n[i] = float(r_line[0])
        ne[i] = float(r_line[1])*pow(10, -20)
        ne_er[i] = float(r_line[2])*pow(10, -20)
        te[i] = float(r_line[3])/1000
        te_er[i] = float(r_line[4])/1000
        i += 1

    profiles['psi_normal'] = psi_n
    profiles['electron_density(10^20/m^3)'] = ne
    profiles['density error(10^20/m^3)'] = ne_er
    profiles['electron_temperature(KeV)'] = te
    profiles['temperature error(10^20/m^3)'] = te_er
    return profiles

def tanh(r,r0,h,d,b,m):
    return b+(h/2)*(np.tanh((r0-r)/d)+1) + m*(r0-r-d)*np.heaviside(r0-r-d, 1)

def expfit(x,A,l):  #Removed vertical displacement variable B; seemed to cause 'overfitting'
    return A*np.exp(l*x)


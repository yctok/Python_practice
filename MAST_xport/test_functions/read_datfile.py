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

# def find_time(mastfile_loc):
#     timeid = mastfile_loc[mastfile_loc.rfind('_')+1:mastfile_loc.rfind('.')]
#     return timeid

# print(find_time('yag_27205_275.dat'))


# def check_tail(mastfile_loc):
#     if mastfile_loc[-4:] == '.dat':
#         return True
#     else:
#         return False

# print(check_tail('yag_27205_275.dat'))



# def check_shotnum(gfile_loc, shotnum = None):
#     if shotnum is None:
#         try:
#             shotnum = int(gfile_loc[gfile_loc.rfind('g')+1:gfile_loc.rfind('.')])
#             return shotnum
#         except:
#             return None
    
# print(check_shotnum('g027205.00275_efitpp'))    


# import numpy as np
# import SOLPSutils as sut

# def calcPsiVals(gfile_loc, shift):
#     """
#     Call b2plot to get the locations of each grid cell in psin space

#     Saves the values to dictionaries in self.data['solpsData']
#     """

#     """
#     Find grid corners first:
#       0: lower left
#       1: lower right
#       2: upper left
#       3: upper right

#     Average location of cells 0 and 2 for middle of 'top' surface, 
#     which is the top looking at outboard midplane
#     Don't average over whole cell, dR << dZ at outboard midplane 
#     and surface has curvature, so psin will be low

#     jxa = poloidal cell index for the outer midplane
#     crx = radial coordinate corner of grid [m]
#     cry = vertical coordinate corner of grid [m]
#     writ = write b2plot.write file
#     f.y = plot against y
#     """
#     g = sut.loadg(gfile_loc)
#     d = float(shift)
    
#     dR = g['rdim'] / (g['nw'] - 1)
#     gR = []
#     for i in range(g['nw']):
#         gR.append(g['rleft'] + i * dR + d)

#     gR = np.array(gR)
#     return gR

# if __name__ == '__main__':
#     radius_cor = calcPsiVals('g027205.00275_efitpp',0)
#     m_cor = calcPsiVals('g027205.00275_efitpp',1)






    
    



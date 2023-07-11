# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:13:25 2023

@author: Yi-Cheng
"""

import latexexpr
v1 = latexexpr.Variable('H_{ello}',3.25,'m')
print('$$ %s $$'%v1)
v2 = latexexpr.Variable('W^{orld}',5.63,'m')
print('$$ %s $$'%v2)
e1 = latexexpr.Expression('E_{xample}',v1+v2,'m')
print('$$ %s $$'%e1)
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 14:19:37 2023

@author: Yi-Cheng
"""

"from SymPy Tutorial (2022) by Mr. P Solver"
"https://www.youtube.com/watch?v=1yBPEPhq54M&list=PLkdGijFCNuVm4IfZlsZPEt4fPJHfl-0g5"

import sympy as smp
import numpy as np
import matplotlib.pyplot as plt

# x = smp.symbols('x')
# print(smp.cos(x))
# a = smp.sin(x)
# print(a- x**2)
# y = x**2 + 4*x + 3
# z = y**2
# # z.factor()
# print(z.factor())
# # z.expand()
# print(z.expand())
# an = smp.solve(z,x)
# print(an)


"See if factor function can factor with imaginary number"

# b1 = x**2 + 1
# b2 = b1**2
# print(b2.factor())
# anb = smp.solve(b2, x)
# print(anb)  
"imaginary solution can be found!"

"specify our symbols"

"define solution as real"
# x = smp.symbols('x', real=True)
# ax = smp.solve(x**2 + 1, x)
# print(ax)

"define solution as real and positive"
# x = smp.symbols('x', real=True, positive=True)
# ax1 = smp.solve(x + 4, x)
# print(ax1)

"Define many variable at once"

x, y, z = smp.symbols('x y z')
F = x**2 + smp.sin(z)*y
print(F)

x_sols = smp.solve(F, x)
print(x_sols)

y_sols = smp.solve(F, y)
print(y_sols)

z_sols = smp.solve(F, z)
print(z_sols)

"assign a number to the symbolic solution"

"For ploting"
expr = z_sols[0]
print(expr)
expr_f = smp.lambdify([x, y], expr)
print(expr_f(1, 2))

"in general"
ex1 = F.subs([(y, 3), (z, smp.pi/2)])
print(ex1)
ex2 = F.subs([(y, smp.cos(z)),(z,y)])
print(ex2)
"first substitute y then substitute z"


"plot the symbolic solution"

x_num = np.linspace(0,1,100)
y_num = 2
plt.plot(x_num, expr_f(x_num, y_num))
plt.show()










 



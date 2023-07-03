# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 06:01:57 2023

@author: user
"""

import sympy as smp
import numpy as np
import matplotlib.pyplot as plt

"""
A falling object encounters a moving platform accelerating upwards:
    1. object: h0(t) = h0 - v0*t -(1/2)*g*t^2
    2. platform: hp(t) = vp*t + (1/2)*q*t^2

find the initial velocity v0 such that when the object and platform collide,
they are moving at the same speed

we need to solve for v0 and t in the two equation
    1. h0(t) - hp(t) = 0
    2. dh0(t)/dt + dhp(t)/dt = 0

"""

"Define expressions"

t, h0, v0, g, vp, q = smp.symbols('t h_0 v_0 g v_p q', real=True)
h0t = h0 - v0*t - smp.Rational(1,2)*g*t**2
dh0dt = -v0 + g*t
hpt = vp*t + smp.Rational(1,2)*q*t**2
dhpdt = vp + q*t


"Define equations"

eq1 = h0t - hpt
eq2 = dh0dt + dhpdt

"solve the equations"

ans = smp.solve([eq1, eq2], [t, v0])
print(ans)
t_collide, v_initfall = smp.solve([eq1, eq2], [t, v0])[1]
print(t_collide)
print(v_initfall)

a = dh0dt.subs([(t,t_collide),(v0,v_initfall)]).simplify()
print(a)
b = dhpdt.subs([(t,t_collide),(v0,v_initfall)]).simplify()
print(b)


h0t_f = smp.lambdify([t, h0, v0, g], h0t)



"plot the symbolic solution"

t_num = np.linspace(0,5,100)
h0_num = 400
g_num = 9.8
q_num = 5
vp_num = 8
v0_num = v_initfall.subs([(vp, vp_num), (g, g_num),(h0, h0_num), (q, q_num)])
print(v0_num)
t_c = t_collide.subs([(vp, vp_num), (g, g_num),(h0, h0_num), (q, q_num)])
print(t_c)
dh0t_c = dh0dt.subs([(t, t_c),(v0, v0_num), (g, g_num)])
print(dh0t_c)
h0t_c = h0t.subs([(t, t_c), (h0, h0_num), (v0, v0_num), (g, g_num)])
# a = h0t_f(1, h0_num, v0_num, g_num)
# print(a)
# plt.scatter(1, h0t_f(1, h0_num, v0_num, g_num))
# plt.scatter(1, t_c)
for tt in t_num:
    if tt <= t_c:
        plt.scatter(tt, h0t_f(tt, h0_num, v0_num, g_num), color='blue')
    elif tt > t_c:
        plt.scatter(tt, h0t_f(tt, h0t_c, dh0t_c, g_num), color='red')
        
plt.show()




# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 21:53:19 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import quad
from scipy.integrate import solve_ivp

"""
Purpose of interpolation

Given x_data[...] and y_data[...], we want to create a function y= f(x),
so we can plug in any value of x and obtain a corresponding value for y

The basic form of interpolation connecting the point in y_data with a straight line,
then y= f(x) = connecting of straight lines
(with different starting points and slopes depends on the y_data)

"""

x_data = np.linspace(0,5,5)
y_data = x_data**2

plt.scatter(x_data, y_data)
plt.show()

"Picture linear interpolation as a plot"

plt.plot(x_data, y_data, 'o--')
plt.show()

"make a interpolation function"
"*y_f is a function, give an x, it will return the interpolate y value"

y_f = interp1d(x_data, y_data, 'linear')
y_new = y_f(1.5)
print(y_new)

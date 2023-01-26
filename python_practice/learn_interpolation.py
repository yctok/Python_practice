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

x_data = np.linspace(0,5,5)
y_data = x_data**2

plt.scatter(x_data, y_data)
plt.show()

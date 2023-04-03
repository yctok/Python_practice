# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:44:55 2023

@author: Yi-Cheng
"""

"Learn SymPy"

"SymPy, Youtube, by TM Quest"

"Introduction to SymPy"

import sympy as sp
import math
import numpy as np
import matplotlib as plt

"A symbolic computation system"

# a = sp.sqrt(2)
# print(a)

# b = sp.sqrt(2)**2
# print(b)

# c = math.sqrt(2)
# print(c)

# d = math.sqrt(2)**2
# print(d)

"Defining Symbols"

"learn how to factor and expand polynomials in SymPy"

"Functions and attributes in this lecture"
"""
sp.Symbol(): Defines a new symbol
sp.symbols(): Defines multiple new symbols
sp.factor(): Factor an expression (for example a polynomial)
sp.expand(): Expand and expression
sp.cos(): The cosine function
sp.sin(): The sine function

"""

"Make a symbol"

x = sp.Symbol('x')
print(x)

y = 2*x + 5
print(y)

z = 2*x + x + 4
print(z)




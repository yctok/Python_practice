# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 11:48:25 2023

@author: user
"""

"""
"Decorators"


"""

def outer_function():
    message = 'Hi'
    
    def inner_function():
        print(message)
    return inner_function()

outer_function()
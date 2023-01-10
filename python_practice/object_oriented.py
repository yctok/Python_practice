# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 21:48:37 2023

@author: user
"""

"""
In this code, we will learn about object-oriented programming

"""

item1 = 'phone'
item1_price = 100
item1_quantity = 5
item1_price_total = item1_price * item1_quantity

print(type(item1))#str
print(type(item1_price))#int
print(type(item1_quantity))#int
print(type(item1_price_total))#int

"""
even though we give all the variables item1 in front of them,
for python, they are four variables with different data type

We want to create a framework that price and quantity are related to our item
"""

"Create a class"


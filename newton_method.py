# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 22:10:13 2022

@author: user
"""

"""
The goal:
    find the sqare root of a number using Newton's method with argparse

"""


#find sqare root of x
#g the initial guess

import argparse

parser = argparse.ArgumentParser(description='solve the square root of a number')
parser.add_argument('-x','--number', type= float, metavar='', required=True, help='number to find square root for')
parser.add_argument('-g','--guess', type= float, metavar='', required=True, help='first guess of the square root')
args = parser.parse_args()



def square_root(number, guess):
    g1 = guess
    g2 = 0
    x = number
    delta = abs(g1 - g2)
    while delta >= 0.04:
        g2 = g1
        g1 = 0.5*(g1 + x/(g1))
        delta = abs(g1 - g2)
        print(g1)
        print(g2)
    s_root = g1
    return s_root

   
if __name__ == '__main__':
    result = square_root(args.number, args.guess)
    print(result)


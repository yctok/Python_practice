# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 13:09:24 2023

@author: Yi-Cheng
"""

"""
This is the example program to practice argparse 
This program calculates the volume of a cylinder with a given height and radius
Cylinder volume = (pi)*(height)*(radius**2)

Double Star or (**) is one of the Arithmetic Operator 
(Like +, -, *, **, /, //, %) in Python Language. 
It is also known as Power Operator.

Arithmetic operators priorities order in Decreasing Mode:
() >> ** >> * >> / >> // >> % >> + >> -

"""
import math
import argparse

parser = argparse.ArgumentParser(description='calculate the volume of a cylinder')
parser.add_argument('-r','--radius', type=int, metavar='', required=True, help='radius of a cylinder')
parser.add_argument('-H','--height', type=int, metavar='', required=True, help='height of a cylinder')
group = parser.add_mutually_exclusive_group()
group.add_argument('-q', '--quiet', action='store_true', help='print quiet')
group.add_argument('-v', '--verbose', action='store_true', help='print verbose')
args = parser.parse_args()

"""
'--radius' is the long notation and '-r' is the shorthand notation
'-h' has been reserve for help
With -r, -h, radius and height becomes optional position and can change order with the flags

In these parse arguments
parser.add_argument(radius', type=int, help='radius of a cylinder')
parser.add_argument('height', type=int, help='height of a cylinder')
radius and height are positional arguments
we have to assign them at the right position to get the answer

we add metavar='' (metavar equals to empty string), 
so the help explaination on the terminal is more clear

If we run the Cylinder_volume function with only one variable, say radius,
we will give us: unsupported operand type(s) for *: 'float' and 'NoneType',
which is not clear for us what is happening
Thus, we add required=True to each argument
When we run the Cylinder_volume function with only one variable,
it will state that we miss a variable,
ex: argparse_practice.py: error: the following arguments are required: -r/--radius

"Options to long or short hand answers"
we add two options in our function, one is quiet and the other is verbose
quiet gives us a short answer of the cylinder volume
verbose gives us a long answer

"""


def Cylinder_volume(radius, height):
    vol = (math.pi)*(height)*(radius**2)
    return vol

if __name__ == '__main__':
    volume = Cylinder_volume(args.radius,args.height)
    if args.quiet:
        print(volume)
    elif args.verbose:
        print('Volume of a cylinder with radius %s and height %s is %s' % (args.radius, args.height, volume))
    else:
        print('Volume of a cylinder is %s' %volume)
    
"""
*Type: "run argparse_practice --radius (radius) --height (height)" to run the code
or *Type: "run argparse_practice -r (radius) -H (height)" to run the code
It does not matter if we mixed up the order
Type: "run argparse_practice -H (height) -r (radius)" will still give you the right result
*Type: "run argparse_practice -h" to find out the position of each argument

"""
    
    
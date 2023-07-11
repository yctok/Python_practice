# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 15:09:08 2023

@author: user
"""

import argparse

parser = argparse.ArgumentParser(description='calculate the exponential result of certain base')
parser.add_argument('-b', '--base', type= float, required=True, help= 'base for calculate exponential result')
parser.add_argument('-e', '--exponent', type= float, required=True, help= 'exponent of the base')
args = parser.parse_args()


def exponential(base, exponent):
    ans = (base)**(exponent)
    return ans

if __name__ == '__main__':
    result = exponential(args.base, args.exponent)
    print(result)
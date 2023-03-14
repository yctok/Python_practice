# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 14:13:14 2023

@author: user
"""

import SOLPSxport_dr as sd
import argparse


parser = argparse.ArgumentParser(description='Generate new b2.transport.inputfile files for SOLPS',
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-g', '--gfileloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
parser.add_argument('-p', '--profilesloc', help='location of profs_*.pkl saved profile file', type=str, default=None)
parser.add_argument('-s', '--shotnum', help='shot number; default = None', type=str, default=None)
parser.add_argument('-t', '--timeid', help='time of profile run; default = None', type=str, default=None)
parser.add_argument('-r', '--runid', help='profile run id; default = None', type=str, default=None)
parser.add_argument('-i', '--tiratiofile', help='File location for Ti/TD ratio; default = None', type=str, default=None)
parser.add_argument('-d', '--tdfileloc', help='File location for TD; default = None', type=str, default=None)
parser.add_argument('-f', '--fractional_change', help='Fractional change to transport coefficients; default = 1',
                    type=float, default=1)
parser.add_argument('-sh', '--shift', help='shift of major radius; default = 1', type=float, default=1)
args = parser.parse_args()

if __name__ == '__main__':
    sd.increment_run(gfile_loc=args.gfileloc, profiles_fileloc=args.profilesloc)
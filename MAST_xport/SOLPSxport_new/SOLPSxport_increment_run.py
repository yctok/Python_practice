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
args = parser.parse_args()

if __name__ == '__main__':
    sd.increment_run(gfile_loc=args.gfileloc, profiles_fileloc=args.profilesloc)
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 18:59:05 2023

@author: user
"""

"""
If we comment out the main() function in name_main_practice
    and with one line code: import name_main_practice, 
    the result is: fist module's name is name_main_practice

If we didn't comment out the main() function in name_main_practice
    and with one line code: import name_main_practice, 
    the result is: print "fist module's name is name_main_practice" twice


"""

import name_main_practice

"""
In name_main_practice_2, it has it's own __name__ recorded,
    the print result is: Second module's name is __main__

"""

print("Second module's name is {}".format(__name__))
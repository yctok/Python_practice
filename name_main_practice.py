# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 14:37:46 2023

@author: Yi-Cheng
"""

# print("fist module's name is {}".format(__name__))


"""
This line: if __name__ == '__main__' can help us check:
    Whether this file is run by python now or is imported to another file

"""

# def main():
#     pass


# if __name__ == '__main__':
#     main()

"Run example: show __name__ == '__main__' "

def main():
    print("fist module's name is {}".format(__name__))


if __name__ == '__main__':
    main()
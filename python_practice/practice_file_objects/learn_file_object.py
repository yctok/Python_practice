# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:35:53 2023

@author: Yi-Cheng
"""

"Open file in python"

"""
Files are recommanded to open with context manager

"""

"Example of opening a file without a context manager"

"""
Because the test.txt file is at the same directory as the python file,
    we only need to give the file name

The defult setting for open function is reading, 'r'
    typical options are: Writing 'w', Appending 'a', Reading+Writing 'r+w'

We need to close a file when we are done using it

name: give the name of the file
mod: give the file's open mod, ex: 'r', 'w', 'r+w'
    
"""

# f = open('test.txt', 'r')
# print(f.name)
# print(f.mode)
# f.close()

"Open the file with context manager"

"""
If we use context manager,i.e. key word 'with', we don't have to write close file,
    instead if the code is executed outside the colon, the file is closed

We still have access to the file object although the file is closed

closed: return boolean value of the statement: This file is closed
read(): read the file, and can only be used in the with colon
        If we print f.read(), it will be all the contents in the file
readlines(): return a list, every element in the list is a line in the file
readline(): return a string, it is a line in the file.
            Every time we execute this function,
                it will return the next upcoming line in the file.
            readline() ends with an empty line by defult, we can use:
                end= '' to change it to an empty string 

"""

# with open('test.txt','r') as f:
    # f_contents_a = f.read()
    # print(f_contents_a)
    # f_contents_b = f.readlines()
    # print(f_contents_b)
    # f_contents_c = f.readline()
    # print(f_contents_c, end= '')
    # f_contents_c = f.readline()
    # print(f_contents_c, end= '')
   
# print(f.closed)

"Example for reading a large file"

# with open('test.txt','r') as f:
#     for line in f:
#         print(line, end='')
    
"Use read with specified characters"

"""
In read() function, we can specify how many character to be read.
    If we copy the same line and paste under it,
        it will continue to read the file and give the rest of it
    If the number of the characters we ran have exceeded the total number of
        characters in the file, read will return an empty string

"""

# with open('test.txt','r') as f:
#     f_contents = f.read(50)
#     print(f_contents, end= '')
#     f_contents = f.read(50)
#     print(f_contents, end= '')

"load large file using read"

"""
Change to end= '*', so we can see how the characters are printed

"""

# with open('test.txt','r') as f:
#     read_character_size = 10
#     f_contents = f.read(read_character_size)
    
#     while len(f_contents) > 0:
#         print(f_contents, end ='*')
#         f_contents = f.read(read_character_size)

"find the current reading position"

"""
tell(): return the current reading position

"""

# with open('test.txt','r') as f:
#     read_character_size = 10
#     f_contents = f.read(read_character_size)
#     print(f.tell())


"Repeat the last time reading"

"""
seek(): use the position assigned to indicate the text reader where to start reading

"""

# with open('test.txt','r') as f:
#     read_character_size = 10
#     f_contents = f.read(read_character_size)
#     print(f_contents, end='')
    
#     f.seek(5)
    
#     f_contents = f.read(read_character_size)
#     print(f_contents)

"Write a file"

"""
If we want to write a file, while we type: open ('test.txt','r'),
    python will tell us that the file is not writable

If we type pass after with open ('test.txt','w') as f:, 
    the test.txt file will still be created

If we type f.write('Test') twice, we will write down 'Test' in our file twice.

Seek in write will replace from the position where we start assign it.

"""


# with open ('test2.txt','w') as f:
#     f.write('Test')
#     f.seek(0)
#     f.write('Play')
    
    
# with open ('replacetest.txt','w') as f:
#     f.write('Test')
#     f.seek(0)
#     f.write('P') 

"""
When we deal with multiple files, 
    a better way to do this is to nest it with one another
    
In the code below, we open our original text file, and create a new file,
    run a for loop to copy all the lines in the original text file to the new file
    

"""

with open('test.txt','r') as rf:
    with open('test_copy.txt','w') as wf:
        for line in rf:
            wf.write(line)
    


    


    
    

    


# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 12:19:17 2023

@author: user
"""

"Tuples"



"""
"Tuples as coordinates"

*Tuples are similar to list.
*A tuple can be created in python using the parentheses (round brackets)
*Tuples are immutable, i.e. once they are created, 
    they cannot be changed or modified,
    i.e. we cannot add or erase any element from it
*If we try to change an element in the tuple,this is the error we will get:
    'tuple' object does not support item assignment
*Tuples are often used for things that does not change, ex: coordinates
*Same with list, tuples index number starts with 0

"""

# coordinates = (4, 5)
# print(coordinates[0])

"list of tuples"

# list_coordinate = [(3, 7), (9, 12), (36, 45)]

"function in python"

"""
A function is a collection of codes that performs certain task
*Use def to tell python that we are going to create a function
*Type open and close parentheses and a colon after the name of the function
*Indent the codes that are used for the function

"""

"""
The function created here is to say hi to the user
"""

# def say_hi(name, age):
#     print("Hello, " + name, "you are " + age)

"""
*The code inside a function is not going to get executed, 
    unless we specifically execute it,i.e. calling a function
*The position of "def" in a code does not matter, 
    python will look for "def" in the order of functions that get executed
*Specify parameter in the function so that we can pass information in it
*Make sure all the element are printed with the same data type

"""

# print("Top")
# say_hi("Mike", "25")
# say_hi("Steve", "35")
# print("Bottom")

"""function accept parameters in all kinds of variables,
    i.e. strings, numbers, booleans

"""


# def say_hi_2(name, age):
#     print("Hello, " + name, "you are " + str(age))

# print("Top")
# say_hi_2("Mike", 25)
# say_hi_2("Steve", 35)
# print("Bottom")

"""
"return statement"
The keyword: return allow functions to send back information
*return also breaks out of the function, 
    so any code type after return in the function will not be executed

"""

"cube a number"

# def cube_0(number):
#     return number*number*number

# print(cube_0(3))

# def cube_1(number):
#     print("Let's have fun")
#     return number**3
    
# print(cube_1(4))

# def cube_2(number):
#     cu = number**3
#     return cu
#     print("Let's have fun")

# print(cube_2(4))
# result = cube_2(5)
# print(result)


"""
"If statement"
*A special struction to let the code make decisions, 
    and allow the code to respond to input data
*In if statement, if certain condition is true, we do something,
    or we go to check other conditions



"""

# is_male = False


# if is_male == True:
#     print("This works")
    
"Short hand version for boolean variable"

# if is_male:
#     print("You are a male")
# else:
#     print("You are not a male")

"add another variable"

# is_tall = True

# "'or' operator"

# if is_male or is_tall:
#     print("You are a male or tall or both")
# else:
#     print("You are neither male nor tall")

"'and' operator"

# if is_male and is_tall:
#     print("You are a tall male")
# else:
#     print("You are either not a male or not tall or neither")

"else if (elif)"

# if is_male and is_tall:
#     print("You are a tall male")
# elif is_male and not(is_tall):
#     print("You are a male but not tall")
# elif not(is_male) and is_tall:
#     print("You are not a male but tall")
# else:
#     print("You are neither male nor tall")


"Combine all"


# is_male = False
# is_tall = True

# if is_male and is_tall:
#     print("You are a tall male")
# elif is_male and not(is_tall):
#     print("You are a male but not tall")
# elif not(is_male) and is_tall:
#     print("You are not a male but tall")
# else:
#     print("You are neither male nor tall")
    

"If statement and comparisons"

"function that find the largest number"
"""
*'>=' is a comparison operator
*statement of comparison returns a boolean (true or false) value

"""

# def max_num(num_1, num_2, num_3):
#     if num_1 - num_2 >= 0 and num_1 - num_3 >= 0:
#         print(str(num_1) + " is the largest")
#     elif num_1 - num_2 >= 0 and not(num_1 - num_3 >= 0):
#         print(str(num_3) + " is the largest")
#     elif not(num_1 - num_2 >= 0) and num_1 - num_3 >= 0:
#         print(str(num_2) + " is the largest")
#     elif not(num_1 - num_2 >= 0) and not(num_1 - num_3 >= 0) and num_2 -num_3 >=0:
#         print(str(num_2) + " is the largest")
#     else:
#         print(str(num_3) + " is the largest")
        
# max_num(5,11,4)

# def max_number(num1, num2, num3):
#     if num1 >= num2 and num1 >= num3:
#         return num1
#     elif num2 >= num1 and num2 >= num3:
#         return num2
#     else:
#         return num3

# print(max_number(41, 25, 37))
# result = max_number(24, 65, 31)
# print(str(result))
# print(result)

# def max_num_2(num_1, num_2, num_3):
#     if num_1 - num_2 >= 0 and num_1 - num_3 >= 0:
#         print(str(num_1) + " is the largest")
#     elif num_2 - num_1 >= 0 and num_2 - num_3 >= 0:
#         print(str(num_2) + " is the largest")
#     else:
#         print(str(num_3) + " is the largest")
        
# max_num_2(26,37,55)

"""
"Comparison operators"
'!=': not equal to
'==': equal to
'<=': less than or equal to
'>=': larger than or equal to
'<': less than
'>': larger than

"""

"Building a calculator"

# import argparse

# parser = argparse.ArgumentParser(description='basic calculator')
# parser.add_argument('-a', '--number_1',type= int, required=True, help='first number for calculation')
# parser.add_argument('-d', '--operator',type= str, required=True, help='define the operator')
# parser.add_argument('-b', '--number_2',type= int, required=True, help='second number for calculation')
# args = parser.parse_args()

# def calculator(number_1, operator, number_2):
#     cal = [number_1, number_2, operator]
#     return cal
    

# if __name__ == '__main__':
#     result = calculator(args.number_1, args.operator, args.number_2)
#     if args.operator == "'+'":
#         print("I see +")
#         val = (result[0])+(result[1])
#         print(val)
#     elif args.operator == "'-'":
#         print("I see -")
#         val = (result[0])-(result[1])
#         print(val)
#     elif args.operator == "'*'":
#         print("I see *")
#         val = (result[0])*(result[1])
#         print(val)
#     elif args.operator == "'/'" and result[1] != 0:
#         print("I see /")
#         val = (result[0])/(result[1])
#         print(val)
#     elif args.operator == "'/'" and result[1] == 0:
#         print("we have trouble dividing zero")
#     elif args.operator == "'**'":
#         print("I see **")
#         val = (result[0])**(result[1])
#         print(val)
#     else:
#         print("Invalid operator")


"Debug tool"
    
    # if args.operator == "'*'":
    #     print(result)
    # else:
    #     print("there is a bug")
    #     print(result)

"Clear version"

# import argparse

# parser = argparse.ArgumentParser(description='basic calculator')
# parser.add_argument('-a', '--number_1',type= float, metavar= '', required=True, help='first number for calculation')
# parser.add_argument('-d', '--operator',type= str, metavar= '', required=True, help='define the operator')
# parser.add_argument('-b', '--number_2',type= float, metavar= '', required=True, help='second number for calculation')
# args = parser.parse_args()

# def calculator(number_1, operator, number_2):
#     cal = [number_1, number_2, operator]
#     return cal
    

# if __name__ == '__main__':
#     result = calculator(args.number_1, args.operator, args.number_2)
#     if args.operator == "'+'":
#         print((result[0])+(result[1]))
#     elif args.operator == "'-'":
#         print((result[0])-(result[1]))
#     elif args.operator == "'*'":
#         print((result[0])*(result[1]))
#     elif args.operator == "'/'" and result[1] != 0:
#         print((result[0])/(result[1]))
#     elif args.operator == "'/'" and result[1] == 0:
#         print("we have trouble dividing zero")
#     elif args.operator == "'**'":
#         print((result[0])**(result[1]))
#     else:
#         print("Invalid operator")

"dictionaries"

"""
Dictionaries allow data to be stored in key-value pairs
key -> words, value -> definition and explanation

*Dictionaries are created with a curly bracket (braces),
    all key-value pairs are stored in it
*key and value have to be a one to one relation
*get(): return value correspond to the assigned key
    **get() returns a value 'None' when the assigned key's value is not found


"""

"Store key-value pair in strings"

# Spring_Conversion = {
#     "Feb": "February",
#     "Mar": "March",
#     "Apr": "April",
#     }

# print(Spring_Conversion["Apr"])
# print(Spring_Conversion.get("Mar"))
# print(Spring_Conversion.get("Jul"))
# print(Spring_Conversion.get("Jul", "Not a valid key"))

"Store key-value pair with numbers and strings"

"Strings as key"

# Spring_Conversion = {
#     "Feb": 2,
#     "Mar": 3,
#     "Apr": 4,
#     }

# print(Spring_Conversion["Apr"])
# print(Spring_Conversion.get("Mar"))
# print(Spring_Conversion.get("Jul"))
# print(Spring_Conversion.get("Jul", "Not a valid key"))


"Strings as value"

# Spring_Conversion = {
#      2: "February",
#      3: "March",
#      4: "April",
#     }

# print(Spring_Conversion[4])
# print(Spring_Conversion.get(3))
# print(Spring_Conversion.get(7))
# print(Spring_Conversion.get(7, "Not a valid key"))

"Mixed key"

Spring_Conversion = {
     2: "February",
     3: "March",
     4: "April",
     "Feb": 2,
     "Mar": 3,
     "Apr": 4,
    }

print(Spring_Conversion[4])
print(Spring_Conversion.get(3))
print(Spring_Conversion["Apr"])
print(Spring_Conversion.get("Mar"))
print(Spring_Conversion.get(7))
print(Spring_Conversion.get(7, "Not a valid key"))

"While loop"









    
    
    


    
    

  









    





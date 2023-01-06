# -*- coding: utf-8 -*-
"""
Spyder Editor

python practice

"get familiar with python"

"""

"""
practice print function

"""

# print("Hello world!!")
# print("  /")
# print(" /")
# print("/")

"""
print with a string in it

"""

# Character_name = "Peter"
# Character_age = "35"
# print("Once, there was a man named " + Character_name + ",")
# print("he is " +Character_age+ " years old")

# Character_name ="Tom"
# print("However, he now prefer to be called by " +Character_name+ ".")

"""

"working with strings in python"
String is way to store words and symbles
we can also store numbers
we can work with True/False value (Boolean Values)

"""

# print("Giraffe\nAcademy")
# print("Giraffe\"..\"Academy")
# print("Giraffe\..\Academy")

# phrase = "william and mary tokamak"
# print(phrase)

"""
"concatenation"
practice to connect a series of strings
"""

# phrase = "william and mary tokamak"
# print(phrase + " is cool")


"""

"operate with functions"
practice using functions associated with strings praperties
lower: return lower case of each character in the strings
upper: return upper case of each character in the strings
islower: return Boolean values of a string, lower -> true
isupper: return Boolean values of a string, upper -> true
    *The isupper result of "PracTice" is false
    *The islower result of "PracTice" is false
len: return the length of a string
    *The length calculate by len function includes space
[0]: the index number assigned to each character in a string, 
    ex: 0 is the first character in a string
    *If we assign a number longer than the length of a string, there will be an error
    *print index number have to start from smaller number to larger number for both positive and negative numbers
        **If print the index number with wrong order, it will show nothing
index: return the number assigned to the input character
    *index returns the first index number of the charactor it assigned to find

replace: return a new string with replaced characters
    *If there several same words in a string, can we replace at once? Yes!
    *' food[3] = "icecream" ' can not be used in string variable


curious:
How to replace a word in a code?
 
"""

phrase = "william and mary tokamak"
# print(phrase + " is cool")
fn = "WE ARE THE CHAMPIONS"
fm = "PracTice"
tn = "P A Q"
# print(fn.lower())
# print(fn.upper())
# print(fn.upper().isupper())
# print(fm.isupper())
# print(fm.islower())
# print(len(tn))
# print(len(phrase))
# print(tn[0])
# print(phrase[-2:-5])
# print(phrase[2:4])
# print(fn[5])
# print(fn.index("M"))
# print(fn.replace("CHAMPIONS","SUPERSTARS"))
# print(fn.replace("E","BB"))

"""

"Working with numbers in python"
"basic math calculation"
%: mod

"""

# print(2)
# print(2.0986)
# print(-2.0986)
# print(2 + 3.7)
# print(4*6)
# print(10/5)
# print(4-3)
# print(4*2+7)
# print(4*(2+7))
# print(45 % 2)

"""
"math calculation 2"
str(): change input float number into string
python can print an array
abs: return the absolute value of a number
pow(a,b): return the value of a to the power of b
max(a,b,...): return the maximum value among a,b,...
min(a,b,...):  return the maximum value among a,b,...
round(number, digits): return a rounded number with accuracy to assigned digits

"""

# my_num = 7
# print(str(my_num) + " is my favorite number.")
# my_ar = [7, 2, 5]
# print(my_ar)
# my_number = -9
# print(abs(my_number))
# print(pow(2,3))
# print(pow(3,-1))
# print(pow(-7,4))
# print(max(-20,5))
# print(min(-20,5))

"use both max, min, and pow"
# print(max(pow(3,4),pow(4,3)))
# print(min(pow(3,4),pow(4,3)))
# print(round(4.2))

"""
"import python math module"
math.floor: round the input number to the closest smaller integer
math.ceil:  round the input number to the closest larger integer
math.sqrt: return square root of the imput number

"""

# import math

# print(math.floor(3.7))
# print(math.ceil(3.7))
# print(math.sqrt(81))

"Small projects"

"1. interaction with user"

# print("Enter your name: ")
# name = input()
# print("Hello, " + name)

"2. build a calculator"

# num1 = 7.2
# num2 = 8.6
# result1 = int(num1) + int(num2)
# result2 = float(num1) + float(num2)
# re = [result1, result2]
# print(re[0])
# print(re)

"""

"lists"

"""


# friends = ["Kevin","Rob","Derek"]
# print(friends[1])
# print(friends[-1])
# print(friends[1:])
# print(friends[-2:])
# print(friends[-3:-2])

# flist = ["Kevin","Rob","Derek", "Oscar", "Pravin"]
# print(flist[1:4])

# food = ["chocolate","cookie","candy","bread"]
# food[3] = "icecream"
# print(food)
# print(food[3])


"list function"

lucky_number = [4, 8, 6, 10, 15]
friends =  ["Kevin","Rob","Derek", "Oscar", "Pravin"]
fr =  ["Kevin","Rob","Derek", "Oscar", "Pravin"]
friends.extend(lucky_number)
print(friends)
friends.extend(fr)
print(friends)











# -*- coding: utf-8 -*-
"""
Spyder Editor

python practice
"""
"get familiar with python"

'''
print("Hello world!!")
print("  /")
print(" /")
print("/")

Character_name = "Peter"
Character_age = "35"
print("Once, there was a man named " + Character_name + ",")
print("he is " +Character_age+ " years old")

Character_name ="Tom"
print("However, he now prefer to be called by " +Character_name+ ".")
'''   

'''
String is way to store words and symbles
we can also store numbers
we can work with True/False value

'''


"working with strings in python"

'''

print("Giraffe\nAcademy")
print("Giraffe\"..\"Academy")
print("Giraffe\..\Academy")

phrase = "william and mary tokamak"
print(phrase)

"concatenation"
phrase = "william and mary tokamak"
print(phrase + " is cool")


"operate with functions"
fn = "WE ARE THE CHAMPIONS"
print(fn.lower())
print(fn.upper())
print(fn.upper().isupper())
print(len(fn))
print(len(phrase))
print(phrase[0])
print(fn[5])
print(fn.index("M"))
print(fn.replace("CHAMPIONS","SUPERSTARS"))

'''

"Working with numbers in python"

"basic math calculation"

'''
print(2)
print(2.0986)
print(-2.0986)
print(2 + 3.7)
print(4*6)
print(10/5)
print(4-3)
print(4*2+7)
print(4*(2+7))
print(45 % 2)

'''

"math calculation 2"
'''

my_num = 7
print(str(my_num) + " is my favorite number.")
my_ar = [7, 2, 5]
print(my_ar)
my_number = -9
print(abs(my_number))
print(pow(2,3))
print(pow(3,-1))
print(pow(-7,4))
print(max(-20,5))
print(min(-20,5))
print(max(pow(3,4),pow(4,3)))
print(min(pow(3,4),pow(4,3)))
print(round(4.2))

'''
"import python math"
'''
import math

print(math.floor(3.7))
print(math.ceil(3.7))
print(math.sqrt(81))

"interaction with user"

print("Enter your name: ")
name = input()
print("Hello, " + name)
'''

"build a calculator"

'''
num1 = 7.2
num2 = 8.6
result1 = int(num1) + int(num2)
result2 = float(num1) + float(num2)
re = [result1, result2]
print(re[0])
print(re)

'''

"mad libs game"

"lists"

'''
friends = ["Kevin","Rob","Derek"]
print(friends[1])
print(friends[-1])
print(friends[1:])
print(friends[-2:])
print(friends[-3:-2])

flist = ["Kevin","Rob","Derek", "Oscar", "Pravin"]
print(flist[1:4])

food = ["chocolate","cookie","candy","bread"]
food[3] = "icecream"
print(food[3])

'''

"list function"

lucky_number = [4, 8, 6, 10, 15]
friends =  ["Kevin","Rob","Derek", "Oscar", "Pravin"]
fr =  ["Kevin","Rob","Derek", "Oscar", "Pravin"]
friends.extend(lucky_number)
print(friends)
friends.extend(fr)
print(friends)











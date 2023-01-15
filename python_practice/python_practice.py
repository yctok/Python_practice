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
# fn = "WE ARE THE CHAMPIONS"
# fm = "PracTice"
# tn = "P A Q"
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
list can be created in python using brackets (square brackets)
    *python can print a list
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

"""
"list function"
*The reason why we need extra same array because friends array have changed after we apply .extend function
I have a question that fr = friends.extend(lucky_number) does not work!

list can be created in python using brackets (square brackets)
    *python can print a list

append(): return a new list with a new item added at the bottom
count(): return the number of how many same items are there in the list
extend(): return a new list with another list added at the bottom
index(): return the index (position) number of the assigned item in the list
    *If we assign a item that is not in the list, ex:"Mike", 
    python will return an error: "Mike" is not in the list
insert(): return a new list with a new item added at the assigned position
    *the flag index and element in the insert() function might be a position flag
pop(): return a new list with the last item in the list get removed
remove(): return a new list with the assigned item get removed
    *If there are multiple same elements in the list,
    the element with the smallest index will be removed
reverse(): return the original list with all elements sorted in reversed order
sort(): return the original list with all string elements sorted in alphabetical order
    *for list with number elements, they will be sorted from small to large
    *sort() function does not support for list with both string elements and numbers


clear(): clear out all the items in the list and return a empty list
copy(): return the copy of a list
    *copy result can be used as a variable
"""

"list with all string elements"

# name = ["Tom","Tom","Tom", "Richard","Richard", "May"]
# friends =  ["Kevin","Rob","Derek", "Oscar", "Pravin"]
# print(friends)
# friends.append("Bob")
# print(friends)
# print(name.count("Tom"))
# print(name.count("Richard"))
# friends.extend(name)
# print(friends)
# print(friends.index("Rob"))
# friends.insert(2, "Richard")
# print(friends)
# friends.pop()
# print(friends)
# friends.remove("Derek")
# print(friends)
# friends.reverse()
# print(friends)
# friends.sort()
# print(friends)

# f2 =  friends.copy()
# print(f2)
# f2.clear()
# print(f2)


"list with all number elements"

# lucky_number = [4, 8, 6, 10, 15]
# number = [7, 61, 42, 7, 23, 31, 7, 96]
# lucky_number.append(41)
# print(lucky_number)
# print(number.count(7))
# lucky_number.extend(number)
# print(lucky_number)
# print(lucky_number.index(10))
# lucky_number.insert(3, 21)
# print(lucky_number)
# lucky_number.pop()
# print(lucky_number)
# lucky_number.remove(7)
# print(lucky_number)
# lucky_number.reverse()
# print(lucky_number)
# lucky_number.sort()
# print(lucky_number)

# n2 = [4, 8, 6, 10, 15]
# print(n2)
# n2.clear()
# print(n2)


"list with both strings and numbers"

# friend_w_luck = ["Kevin", 7,"Rob", 4, "Derek", 7, "Oscar", "Pravin", 11, 2, 7, 35]
# fun = ["python", "all", 41]
# friend_w_luck.append(25)
# print(friend_w_luck)
# friend_w_luck.append("May")
# print(friend_w_luck)
# l7 = friend_w_luck.count(7)
# print(l7)
# print(friend_w_luck.count("Rob"))
# friend_w_luck.extend(fun)
# print(friend_w_luck)
# print(friend_w_luck.index(4))
# friend_w_luck.insert(4, "night")
# print(friend_w_luck)
# friend_w_luck.pop()
# print(friend_w_luck)
# friend_w_luck.remove("all")
# print(friend_w_luck)
# friend_w_luck.reverse()
# print(friend_w_luck)

# fwl2 = ["Kevin", 7,"Rob", 4, "Derek", 7, "Oscar", "Pravin", 11, 2, 7, 35]
# print(fwl2)
# fwl2.clear()
# print(fwl2)











# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 12:05:49 2023

@author: user
"""


"Object-Orianted tutorials"

"from Youtube Corey Schafer, https://www.youtube.com/playlist?list=PL-osiE80TeTsqhIuOqKhwlXsIBIdSeYtc"

"Tutorial 1: Class and Instances"

"""
why we use class?
Class allow us to logically group our data and function in way which is easy to reuse,
    and also easy to build upon if needed

For data and functions that are built associated with a class,
    we called them attributes and methods.
    Data -> attributes
    Functions -> methods

A class is framework. In our following example, 
    Employee class define what we need to know about an employee and what action can be done
    An employee is an instance of a class
    Employee() is a class
    emp_1 = Employee() is an instance of a class

"""

"Create an employee class for a company"

"Inefficient way of using class"

# class Employee:
#     pass

# emp_1 = Employee()
# emp_2 = Employee()

# print(emp_1)
# print(emp_2)

# emp_1.first = 'Corey'
# emp_1.last = 'Schafer'
# emp_1.email = 'Corey_Schafer@w&m.com'
# emp_1.pay = 50000

# emp_2.first = 'Yi-Cheng'
# emp_2.last = 'Chuang'
# emp_2.email = 'yichengchuang@w&m.com'
# emp_2.pay = 40000

# print(emp_1.email)
# print(emp_2.email)

"""
"__init__ method"

Every method in python take a instance of the class as their first parameter,
    it is called self in convension

We can also define it as self.first_name = first, 
    it just that keeping them the same is easier

Defining self.first = first is the same as defining emp_1.first = 'Corey',
    because we use __init__  method, 
    python will construct the according relation for us follow the given pattern

The self.email does not have independent variable,
    but it is an independent attribute (data)

"fullname function"

We create a function that can generate a employee's full name

Because there is no new information needed for fullname() function,
    the variable in fullname() function is self only

When we put: '{} {}'.format(emp_1.first, emp_1.last) in the function,
    we need to change emp_1 to the general instance representation, self.

Remember to add () ,i.e. parentheses, at the end of fullname() function,
    because fullname() is a method (function) not a attribute (data)
    If we forget the parentheses, python will tell us its property, a method.

we run fullname() function from class: Employee.fullname(emp_1) and reach the same result,
    which shows that a instance (self) is a variable pass in the fullname() function.

"""

# class Employee:
    
#     def __init__(self, first, last, pay):
#         self.first = first
#         self.last = last
#         self.pay = pay
#         self.email = first + '_' + last + '@w&m.com'
    
#     def fullname(self):
#         return '{} {}'.format(self.first, self.last)



# emp_1 = Employee('Corey', 'Schafer', 50000)
# emp_2 = Employee('Yi-Cheng', 'Chuang', 40000)

# print(emp_1.email)
# print(emp_2.email)

# print('{} {}'.format(emp_1.first, emp_1.last))
# print(emp_1.fullname())
# print(Employee.fullname(emp_1))
# print(emp_2.fullname())

"Tutorial 2: Class variables"

"""
"Class variables"

Class variables are those we can apply for every instance in the class.

apply_raise() function renews all the pay attributes

"raise_rate print example"

Class variables can be accessed through the class as well as instances.

"raise_rate variable examine"

There is no raise_rate store in instances, 
    but python will also search for this variable in the class

"define class variable outside the class"

we can define the class variable outside the class, 
    ex: Employee.raise_rate = 1.04

we can define a instance variable outside the class,
    ex: emp_1.raise_rate = 1.03
    
If we define an instance variable, it will be recorded in the instance's __dict__

It is better to define self.raise_rate in apply_raise,
    (better than define with the class: Employee.raise_rate)
    because we then have freedom to define a class variable or an instance variable
    (also easier for subclass to rewrite it)

"Define number of employee"

Because number of employee is collective information, 
    we define it as a class variable

__init__ method can help us to count the number of employee

"""

# class Employee:
#     num_of_emps = 0
#     raise_rate = 1.04
    
#     def __init__(self, first, last, pay):
#         self.first = first
#         self.last = last
#         self.pay = pay
#         self.email = first + '_' + last + '@w&m.com'
        
#         Employee.num_of_emps += 1
    
#     def fullname(self):
#         return '{} {}'.format(self.first, self.last)
    
#     def apply_raise(self):
#         self.pay = int(self.pay * self.raise_rate)

    
# print(Employee.num_of_emps)

# emp_1 = Employee('Corey', 'Schafer', 50000)
# emp_2 = Employee('Yi-Cheng', 'Chuang', 40000)

# print(emp_1.pay)
# emp_1.apply_raise()
# print(emp_1.pay)

"raise_rate print example"

# print(Employee.raise_rate)
# print(emp_1.raise_rate)
# print(emp_2.raise_rate)

"raise_rate variable examine"

# print(emp_1.__dict__) #no raise_rate variable
# print(Employee.__dict__)

"Define class variable outside the class"

# Employee.raise_rate = 1.05

# print(Employee.raise_rate)
# print(emp_1.raise_rate)
# print(emp_2.raise_rate)

"Define class variable for an instance"

# emp_1.raise_rate = 1.03

# print(emp_1.__dict__)

# print(Employee.raise_rate)
# print(emp_1.raise_rate)
# print(emp_2.raise_rate)

"Number of employee"

# print(Employee.num_of_emps)


"Tutorial 3: Static method and Class method"

"""
The method we written above is a regular method,
    it takes self as its first variable

A class method takes class as its first variable

use class method to set raise_amount = 1.05 is equal to set it like:
    Employee.raise_amt = 1.05
    We can also use a instance: emp_1.set_raise_amt =1.05 (not recommend)

"Class method as an ulternative constructor"

Use class method to provide multiple ways to create our own objects

"Static method"

Static method do not pass in the first variable automatically, so
    they act like regular function but have logical connection with the class

If we do not pass in an instance or a class as a variable in a method we written
    in the class, then that method is a static method


    
"""

# class Employee:
#     num_of_emps = 0
#     raise_amount = 1.04
    
#     def __init__(self, first, last, pay):
#         self.first = first
#         self.last = last
#         self.pay = pay
#         self.email = first + '_' + last + '@w&m.com'
        
#         Employee.num_of_emps += 1
    
#     def fullname(self):
#         return '{} {}'.format(self.first, self.last)
    
#     def apply_raise(self):
#         self.pay = int(self.pay * self.raise_amt)
    
#     @classmethod
#     def set_raise_amt(cls, amount):
#         cls.raise_amt = amount
    
#     @classmethod
#     def from_string(cls, emp_str):
#         first, last, pay = emp_str.split('-')
#         return cls(first, last, pay)
    
#     @staticmethod
#     def is_workday(day):
#         if day.weekday() == 5 or day.weekday() == 6:
#             return False
#         return True
    
    
# emp_1 = Employee('Corey', 'Schafer', 50000)
# emp_2 = Employee('Yi-Cheng', 'Chuang', 40000)


"equivalent code"

# Employee.set_raise_amt(1.05) 
# Employee.raise_amt = 1.05 

"""
Example problem:
    Suppose the employee data is passed in first-last-pay format, instead of
    using split function all the time, we can create a class method to do so.

"""

# emp_str_1 = 'John-Doe-70000'
# emp_str_2 = 'Steve-Smith-30000'


"Inefficient way"

# first, last, pay = emp_str_1.split('-')

# new_emp_1 = Employee(first, last, pay)
# print(new_emp_1.email)
# print(new_emp_1.pay)

"Work with class method"

# new_emp_2 = Employee.from_string(emp_str_2)
# print(new_emp_2.email)
# print(new_emp_2.pay)

"Static methods"

# import datetime
# my_date = datetime.date(2023, 1, 12)

# print(Employee.is_workday(my_date))

"Tutorial 4: Inheritance - Creating Subclasses"

"""
The reason to create subclasses is to inherit functions from its parent class

Put the parent class inherited from in the parenthesis after the name of the subclass



"""

class Employee:
    
    raise_amt = 1.04
    
    def __init__(self, first, last, pay):
        self.first = first
        self.last = last
        self.pay = pay
        self.email = first + '_' + last + '@w&m.com'
        
    
    def fullname(self):
        return '{} {}'.format(self.first, self.last)
    
    def apply_raise(self):
        self.pay = int(self.pay * self.raise_amt)
    
"""
Introduce new attribute to the subclass, create a init method for suclass,
use super() or parent class name to put in all the attribute from parent class

"""     
    
class Developer(Employee):
    raise_amt = 1.10
    
    def __init__(self, first, last, pay, prog_lang):
        super().__init__(first, last, pay)
        # Employee.__init__(self, first, last, pay)
        self.prog_lang = prog_lang
        
"""
Introduce a new class: Manager
Do not try to pass in a data type, like empty list [] in the init method


"""

class Manager(Employee):
    
    def __init__(self, first, last, pay, employees=None):
        super().__init__(first, last, pay)
        # Employee.__init__(self, first, last, pay)
        if employees is None:
            self.employees = []
        else:
            self.employees = employees
    
    def add_emp(self, emp):
        if emp not in self.employees:
            self.employees.append(emp)
    
    def remove_emp(self, emp):
        if emp in self.employees:
            self.employees.remove(emp)
    
    


        
dev_1 = Employee('Corey', 'Schafer', 50000)
dev_2 = Employee('Yi-Cheng', 'Chuang', 40000)

# print(dev_1.email)
# print(dev_2.email)

dev_3 = Developer('Joe', 'Biden', 50000, 'python')
dev_4 = Developer('Rob', 'Behary', 40000, 'Java')

# print(dev_3.email)
# print(dev_4.email)

"""
help function generater information that help user to understand 
what the subclass inherit

"""

# print(help(Developer))

"""
We can change some attributes in subclass and without affecting the parent class

"""

# print(dev_3.pay)
# dev_3.apply_raise()
# print(dev_3.pay)

# print(dev_1.pay)
# dev_1.apply_raise()
# print(dev_1.pay)

" test subclass new attribute"

print(dev_3.email)
print(dev_3.prog_lang)

















        



    
    
        

# Деление отрезка пополам

from sympy import *

a = float(input("Left point is "))
b = float(input("Right point is "))
tolerance = float(input("Tolerance is "))

if a >= b:
    print("Check input data")
    exit()

x = symbols('x')
exprs = eval(input())


while b - a > tolerance:

    xi = a + (b-a)/2
    if sign(exprs.subs(x,a).n()) != sign(exprs.subs(x,xi).n()):
        b = xi
    else:
        a = xi

print(" Root is ", xi,'\n', "Tolerance is", tolerance)

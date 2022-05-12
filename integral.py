import math
from numpy import *

def func(x):
    return math.sin(100*x)*math.exp(-x**2)*math.cos(2*x)

def Simpson_method(f,x0,xf,N):
    funcs = []
    for i in range(2*N+1):
        funcs.append(f(x0+(xf-x0)*i/(2*N)))
    res = funcs[0] - funcs[2*N]
    for j in range(1,N+1):
        res += 4*funcs[2*j-1] + 2*funcs[2*j]
    return (xf-x0)/(6*N)*res

result = Simpson_method(func,0,3,10000)
print(result)
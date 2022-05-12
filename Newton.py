import numpy as np
import sympy as sym

def Newton(x, y):
   if len(x) == 1:
        return d[x[0]]
   return (Newton(x[1:], y) - Newton(x[: len(x) - 1], y)) / (x[len(x) - 1] - x[0])

x = np.random.randint(10, size=3)
y = np.random.randint(10, size=3)

X = sym.symbols('X')

d = dict(zip(x,y))
P = 0
n = 1

for i in range(len(x)):
    P = P + n * Newton(x[0:i + 1], y)
    n = n * (X - x[i])
print(sym.simplify(P))

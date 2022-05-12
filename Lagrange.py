import numpy as np
import sympy as sym

def Lagrange(x,y):
    P = 0
    for j in range(len(x)):
       l = 1
       for i in range(len(x)):
          if i != j:
            l = l * (X - x[i]) / (x[j] - x[i])
       P = P + l * y[j]
    return P

x = np.random.randint(10, size=3)
y = np.random.randint(10, size=3)

X = sym.symbols('X')

p = Lagrange(x,y)

print(sym.simplify(p))

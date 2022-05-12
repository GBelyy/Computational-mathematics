#Метод Лагранжа
import numpy as np
import sympy as sp
import re
"""""
# точки записывать в виде: x1,y1;x2,y2;...;xn,yn
data = re.split(r'[,/;]', input())
X = np.array(data[::2],dtype=float)
Y = np.array(data[1::2],dtype=float)

x = sp.symbols('x')
phi = []
for i in range(len(X)):
    phi_i = 1
    for xi in X:
        if xi == X[i]:
            continue
        expr = (x - xi)/(X[i] - xi)
        phi_i *= expr
    phi.append(phi_i * Y[i])

print(sp.simplify(sum(phi)))
"""""
#Метод Ньютона

data = re.split(r'[,/;]', input())
X = np.array(data[::2],dtype=float)
Y = np.array(data[1::2],dtype=float)

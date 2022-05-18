import numpy as np
from matplotlib import pyplot as plt

#Инициализация
p0 = 100
pn0 = 150
pnL = 50
length = 500
L = 100
tau = 60*60*24

k = 1e-14
m = 1e-3
f = 0.2
cf = 1e-4
r0 = 1e3
pp = 120

r = lambda p: r0*(1 + cf*(p - pp))
t = 0
T = 8640000000


def solveTriangleSlae(a,b,c,f):
    a, b, c, f = tuple(map(lambda k_list: list(map(float, k_list)), (a, b, c, f)))

    alpha = [-b[0] / c[0]]
    beta = [f[0] / c[0]]
    n = len(f)
    x = [0]*n

    for i in range(1, n):
        alpha.append(-b[i]/(a[i]*alpha[i-1] + c[i]))
        beta.append((f[i] - a[i]*beta[i-1])/(a[i]*alpha[i-1] + c[i]))

    x[n-1] = beta[n - 1]

    for i in range(n-1, 0, -1):
        x[i - 1] = alpha[i - 1]*x[i] + beta[i - 1]

    return x

h = length/L
p = np.full(L, p0)
p[0] = pn0
p[-1] = pnL

while t < T:

    a = np.zeros(L)
    b = np.zeros(L)
    c = np.zeros(L)
    d = np.zeros(L)
    for i in range(1, L - 1):
        r_plus = r(p[i]) if p[i] >= p[i + 1] else r(p[i+1])
        r_minus = r(p[i-1]) if p[i-1] >= p[i] else r(p[i])
        c[i] = k*r_plus/(m*h**2)
        b[i] = k*r_minus/(m*h**2)
        a[i] = (-c[i] - b[i] - f*cf*r0/tau)
        d[i] = -f*cf*r0/tau*p[i]
    a[0] = 1
    b[0] = 0
    c[0] = 0
    d[0] = pn0
    a[-1] = 1
    c[-1] = 0
    b[-1] = 0
    d[-1] = pnL
    p = solveTriangleSlae(c, b, a, d) #считем давление через кфты
    t += tau

plt.plot(np.linspace(0, L, L), p)
plt.xlabel('Координата')
plt.ylabel('Время')
plt.title('Зависимость давления от координаты')
plt.show()
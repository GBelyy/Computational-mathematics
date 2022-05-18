import numpy as np
from math import *

lam = 1e-4


def an_sol(x, y, t):
    return cos(pi * x) * sin(5 * pi * y) * exp(-50 * pi * pi * lam * t)


def first_norm_matrix(matrix):
    rows, columns = matrix.shape
    return max([sum([abs(matrix[i][j]) for i in range(rows)]) for j in range(columns)])


def first_norm_vector(vector):
    return max([abs(vector[i]) for i in range(len(vector))])

# для 3 диаг матриц
def method_Zeidel(L, D, U, f, x_0):
    eps = 1e-8
    x = x_0
    A = L + D + U

    while first_norm_vector(np.dot(A, x) - f) > eps:
        right_side = f - np.dot(U, x)
        for i in range(len(x)):
            x[i] = (right_side[i] - np.dot(L[i, :], x)) / D[i][i]
    return x


# Используем метод прогонки
def solveTriagonalSlae(a, b, c, d):
    N = len(a) - 1

    # Вычисляем p1 и q1
    p = [0] + [-c[0] / b[0]]
    q = [0] + [d[0] / b[0]]

    for i in range(1, N):
        znam = a[i] * p[i] + b[i]
        p.append(-c[i] / znam)
        q.append((d[i] - a[i] * q[i]) / znam)

    # Обратный ход:
    yN = (d[N] - a[N] * q[N]) / (a[N] * p[N] + b[N])
    result = [0 for _ in range(N)] + [yN]
    for i in range(N - 1, -1, -1):
        result[i] = p[i + 1] * result[i + 1] + q[i + 1]
    return result


def psi_1(y, t):
    return sin(5 * pi * y) * exp(-50 * pi * pi * lam * t)


def psi_2(y, t):
    return -sin(5 * pi * y) * exp(-50 * pi * pi * lam * t)


def solve_numerical(M):
    # M - число отрезков по икс
    # K - чилос отрезков по игрек
    # N - число отрезков по t
    K = M
    N = 10000

    h = 1 / M  # h_x = h_y
    tau = 1 / N

    # начальное условие
    phi0 = np.zeros((M + 1) * (K + 1), float)
    for m in range(M + 1):
        for k in range(K + 1):
            phi0[m + (M + 1) * k] = cos(pi * h * m) * sin(5 * pi * h * k)
    phi_old = phi0

    A = 25 * lam / h / h
    B = -(52 * lam / h / h + 1 / tau)
    C = lam / h / h

    a = [0] + [A for i in range(1, M)] + [0]
    c = [0] + [A for i in range(1, M)]
    a2 = [0] + [C for i in range(1, M)] + [0]
    c2 = [0] + [C for i in range(1, M)]
    b = [1] + [B for i in range(1, M)] + [1]

    for n in range(1, N + 1):
        if n % 2 == 1:
            phi = [0. for _ in range(M + 1)]

            for j in range(1, K):
                d = [psi_1(h * j, n * tau)] + [-(
                        phi_old[m + (M + 1) * j] / tau + C * phi[m + (M + 1) * (j - 1)] + C * phi_old[
                    m + (M + 1) * (j + 1)]) for m in range(1, M)] + [psi_2(h * j, n * tau)]
                phi += solveTriagonalSlae(a, b, c, d)

            phi += [0. for _ in range(M + 1)]
        else:
            phi = [0 for _ in range((M + 1) * (K + 1))]
            for j in range(K + 1):
                phi[0 + (M + 1) * j] = psi_1(h * j, n * tau)

            for m in range(1, M):
                d = [0] + [
                    -(phi_old[m + (M + 1) * j] / tau + A * phi[m - 1 + (M + 1) * j] + A * phi_old[m + 1 + (M + 1) * j])
                    for j in range(1, K)] + [0]
                result = solveTriagonalSlae(a2, b, c2, d)

                for j in range(K + 1):
                    phi[m+(M + 1) * j] = result[j]

            for j in range(K + 1):
                phi[M + (M + 1) * j] = psi_2(h * j, n * tau)

        phi_old = phi
    return phi_old


def solve_analytical(M):
    K = M
    h = 1 / M
    # аналитическое решение
    phi_an = np.zeros((M + 1) * (K + 1), float)
    for m in range(M + 1):
        for k in range(K + 1):
            phi_an[m + (M + 1) * k] = an_sol(h * m, h * k, 1)
    return phi_an


def diff_anal_and_sol(M):
    return first_norm_vector(solve_numerical(M) - solve_analytical(M))


M = [4, 8, 16, 32, 64, 128, 256, 512]

result = []
for m in M:
    print(m)
    result.append(diff_anal_and_sol(m))
    print(result[-1])
with open("task4result.txt", "w") as file:
    for i in range(len(result)):
        file.write(str(M[i]) + "," + str(result[i]) + "\n")
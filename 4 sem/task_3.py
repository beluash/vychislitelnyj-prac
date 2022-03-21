import matplotlib.pyplot as plt
import numpy as np
import math
import pylab

#заполняем первоначальные данные
with open('input.txt') as f:
    n = int(f.readline())
    b = list(f.readline().split(' '))
    b = [float(b[i]) for i in range(n)]
    eps = float(f.readline())
    x_last = [0.0] * n #начальное приближение
    x_next = [0.0] * n #текущее приближение

with open('matrix.txt') as f1:
    A = [list(map(float, row.split())) for row in f1.readlines()]

#метод релаксации
def RelaxIter(x_last, x_next, eps, A, b, w):
    while ((Nevazka(A, x_next, b) >= eps)):
        for i in range(n):
            s1 = sum(A[i][j] * x_next[j] for j in range(i))
            s2 = sum(A[i][j] * x_last[j] for j in range(i + 1, n))
            x_next[i] = (1 - w) * x_last[i] + w * (b[i] - s1 - s2) / A[i][i]
        #print('x_next: ', x_next)
        x_last = x_next
    return x_next

def Nevazka(A, x, b):
    result = [0.0] * n
    for i in range(n):
        result[i] = abs(b[i] - sum(A[i][j] * x[j] for j in range(n)))
    result = np.sqrt(sum(result[i] ** 2 for i in range(n)))
    #print('nevazka: ', result)
    return result

#метод Якоби
def Jacobi(A, J0, eps):
    max, n1, n2 = FindMax(A)
    while max > eps:
        J0 = JacobiRotation(A, J0, n1, n2)
        max, n1, n2 = FindMax(A)
    lam = [A[i][i] for i in range(n)]
    return lam, J0

def FindMax(A):
    max = abs(A[0][1])
    n1, n2 = 0, 1
    for i in range(n):
        for j in range(n):
            if (abs(A[i][j]) > max) & (i != j): 
                max = abs(A[i][j])
                n1, n2 = i, j
    return max, n1, n2

def JacobiRotation(A, J0, n1, n2):
    J = [[0.0] * n for i in range(n)]
    for i in range(n):
        J[i][i] = 1

    if A[n1][n1] == A[n2][n2]:
        J[n1][n1], J[n1][n2], J[n2][n1], J[n2][n2] = math.cos(math.pi/4), -math.sin(math.pi/4), math.sin(math.pi/4), math.cos(math.pi/4)
    else:
        theta = (A[n1][n1] - A[n2][n2]) / (2 * A[n1][n2])
        if theta != 0: t = (theta/abs(theta))/(abs(theta) + np.sqrt(1 + theta ** 2))
        else: t = 0
        c = 1/np.sqrt(1 + t ** 2)
        s = c * t
        J[n1][n1], J[n1][n2], J[n2][n1], J[n2][n2] = c, -s, s, c
    
    JT = [[J[j][i] for j in range(len(J))] for i in range(len(J[0]))]
    JA = [[sum(a*b for a,b in zip(row, col)) for col in zip(*A)] for row in J]
    JAJT = [[sum(a*b for a,b in zip(row, col)) for col in zip(*JT)] for row in JA]
    J0 = [[sum(a*b for a,b in zip(row, col)) for col in zip(*J0)] for row in J]   
    
    for i in range(n):
        for j in range(n):
            A[i][j] = JAJT[i][j]
    return J0
    
w = 1
J0 = [[0.0] * n for i in range(n)]
for i in range(n):
    J0[i][i] = 1

#print('Приближенное решение системы: ', RelaxIter(x_last, x_next, eps, A, b, w))
result, J0 = Jacobi(A, J0, eps)
print('Собственные числа методом Якоби: ', result)
print('Собственные векторы:', J0)

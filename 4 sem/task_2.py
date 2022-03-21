import matplotlib.pyplot as plt
import numpy as np
import math
import pylab
import random

def Func(x, c):
    result = 0
    for i in range(k+1):
        result  += c[i] * (x ** i)
    return result

#считывание данных
f = open('input.txt', 'r')
a,b = map(float, f.readline().split(' '))
n = int(f.readline())
x, y = [0.0] * (n), [0.0] * (n)
for i in range(n):
    x[i], y[i] = map(float, f.readline().split(' '))
k = int(f.readline())
f.close()

yrandom = [0.0] * (n)
e = 0.2
for i in range(n):
    yrandom[i] =  y[i] + ((-1) ** i) * e * random.uniform(0, n)

c = [0.0] * (k+1) #массив коэффициентов с_j
c1 = [0.0] * (k+1)
b = [0.0] * (k+1) #массив свободных членов справа
b1 = [0.0] * (k+1)
sums = [[0.0] * (k+1) for i in range(k+1)] #суммы при неизвестных коэффициентах
sums1 = [[0.0] * (k+1) for i in range(k+1)]

#---Gauss
def Gauss(x, y, sums, c, b):
    for i in range(k+1):
        for j in range(k+1):
            for l in range(n):
                sums[i][j] += x[l] ** (i+j) 

    for i in range(k+1):
        for j in range(n):
            b[i] += (x[j] ** i) * y[j]

    for i in range(k+1):
        if sums[i][i] == 0:
            for j in range(k+1):
                if (i == j): continue
                if (sums[j][i] != 0 and sums[i][j] != 0):
                    b[j], b[i] = b[i], b[j]
                    for l in range(k+1):
                        sums[j][l], sums[i][l] = sums[i][l], sums[j][l]

    #находим коэффициенты
    for i in range(k+1):
        for j in range(i+1, k+1):
            M = sums[j][i] / sums[i][i]
            for l in range(i, k+1):
                sums[j][l] -= M * sums[i][l]
            b[j] -= M*b[i]
    for i in range(k, -1, -1):
        s = 0
        for j in range(i, k+1):
            s += sums[i][j] * c[j] 
        c[i] = (b[i] - s) / sums[i][i]
    print(c)


Gauss(x, y, sums, c, b)
Gauss(x, yrandom, sums1, c1, b1)

#найти значение в любой точке
while (True):
    point = input("Введите точку: ")
    if (point == "."): 
        break
    print(Func(float(point), c))

#график
xlist = np.arange(x[0], x[n-1], 0.01)
ylistrandom = [Func(v, c1) for v in xlist]
ylist = [Func(v, c) for v in xlist]
plt.scatter(x, y, color="black") #точки
pylab.plot(xlist, ylist , color="blue")
pylab.plot(xlist, ylistrandom , color="red")
pylab.show()

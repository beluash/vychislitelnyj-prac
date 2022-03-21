import matplotlib.pyplot as plt
import numpy as np
import math
import pylab

def Spline(x, A, X):
    n = len(X)
    for i in range(0, n):
        if ((x >= X[i]) & (x <= X[i+1])):
            return A[0][i] + A[1][i]*(x - X[i]) + A[2][i]*(x - X[i])**2 + A[3][i]*(x - X[i])**3
    return 0

#метод трехточечной прогонки
def SetMatrix(A, hx, hy, B, x0, xn): 
    n = len(hx)
    for i in range(0, n-1):
        if (i > 0):
            A[i][i-1] = hx[i] #верхняя диагональ 
        if (i < n-2):
            A[i][i+1] = hx[i+1] #нижняя диагональ
        if (i == 0):
            A[i][i] = 2*(hx[0] + hx[1]) #главная диагональ
            B[i] = 3*(hy[1]/hx[1] - hy[0]/hx[0]) - x0*hx[0]/2
        else:
            if (i == n-2):
                A[i][i] = 2*hx[i] + 1.5*hx[i+1]
                B[i] = 3*(hy[i+1]/hx[i+1] - hy[i]/hx[i]) - 1.5*(xn - hy[i+1]/hx[i+1])
            else:
                A[i][i] = 2*(hx[i] + hx[i+1])
                B[i] = 3*(hy[i+1]/hx[i+1] - hy[i]/hx[i])

def Progonka(A, f):
    n = len(f)
    x, a, b = [0.0] * n, [0.0] * n, [0.0] * n
    a[1] = -A[0][1] / A[0][0]
    b[1] = f[0] / A[0][0]
    for i in range(1, n-1):
        a[i+1] = -A[i][i+1] / (A[i][i-1] * a[i] + A[i][i])
        b[i+1] = (f[i] - A[i][i-1] * b[i]) / (A[i][i-1] * a[i] + A[i][i])
    x[n-1] = (f[n-1] - A[n-1][n-2] * b[n-1]) / (A[n-1][n-1] + A[n-1][n-2] * a[n-1])
    for i in range(n-2, -1, -1):
        x[i] = a[i+1] * x[i+1] + b[i+1]
    return x

#считывание данных
f = open('input.txt', 'r')
a,b = map(float, f.readline().split(' '))
n = int(f.readline())
x, y = [0.0] * (n+1), [0.0] * (n+1)
for i in range(n+1):
   x[i], y[i] = map(float, f.readline().split(' '))
hx = np.diff(x) #шаг сетки
hy = np.diff(y)
x_0 = float(f.readline())
x_n = float(f.readline())
f.close()

A = [[0.0]* (n-1) for i in range(n-1)] #матрица для прогонки
B = [0.0] * (n-1) #правая часть для прогонки
k = [[0.0] * (n+1) for i in range(4)]  #матрица коэффициентов

SetMatrix(A, hx, hy, B, x_0, x_n) #задание матрицы и вектора для прогонки

#заполнение коэффициентов
for i in range(1, n): #заполняем С
    k[2][i] = Progonka(A, B)[i-1]
k[2][0] = x_0/2
k[2][n] = 1.5 * (x_n/hx[n-1] - (hy[n-1]/hx[n-1])/hx[n-1] - k[2][n-1]/3)

for i in range(0, n): #заполняем A, D, B
    k[0][i] = y[i]
    k[3][i] = (k[2][i+1] - k[2][i])/(3*hx[i])
    k[1][i] = hy[i]/hx[i] - hx[i] * (k[2][i+1] + 2*k[2][i]) / 3

#для проверки
def func(x):
    return np.exp(x)

while (True): 
    point = input("Введите точку: ")
    if (point == "."): 
        break
    print(Spline(float(point), k, x))

#график
xlist = np.arange(x[0], x[n], 0.01)
ySlist = [Spline(v, k, x) for v in xlist]
plt.scatter(x, y, color="black") #точки
pylab.plot(xlist, ySlist , color="blue") #синий сплайн
#yFlist = [func(x) for x in xlist]
#pylab.plot(xlist, yFlist, "--", color="red") #красная функция
pylab.show()

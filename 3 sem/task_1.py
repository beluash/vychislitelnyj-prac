import math
import numpy as np
import matplotlib.pyplot as plt

#Наша функция
def func(x):
    return float("%.10f" % (x ** 10))

#Чебышевские узлы
def Chebushev (a, b, j, n):
    return 0.5*(b + a) + 0.5*(b - a) * math.cos(((2*j+1)*math.pi)/(2*(n+1)))

#Полином Лагранажа
def lagrange(x, n, xArray, yArray):
    polynom = 0
    for i in range(n):
        polynom += yArray[i] * coefLagrange(x, i, xArray)
    return polynom

def coefLagrange(x, i, xArray):
    coef = 1
    for j in range(n):
        if i != j:
            coef *= (x - xArray[j]) / (xArray[i] - xArray[j])
    return coef

#Полином Ньютона
def divedDifference(i, j):
    if (i == 0): 
        return yArray[i]
    if (i == 1): 
        return (func(xArray[j + 1]) - func(xArray[j])) / (xArray[j + 1] - xArray[j])
    return (divedDifference(i - 1, j + 1) - divedDifference(i - 1, j)) / (xArray[j + i] - xArray[j])    

def coefNewton(x, i, xArray):
    coef = 1
    for k in range(i):
        coef *= (x - xArray[k])
    return coef

def newton(x, n, xArray, yArray):
    polynom = 0
    for j in range(n):
        polynom += divedDifference(j, 0) * coefNewton(x, j, xArray)
    return polynom

#Погрешность
def sin(x, n, xArray, xArrayGraph, freq):
    return maximum(xArrayGraph, n, freq) / float(math.factorial(n + 1)) * abs(omega(x, n, xArray))

def omega(x, n, xArray):
    w = 1
    for i in range(n):
        w *= (x - xArray[i])
    return w

def derivative (x, order):
    if (order == 0):
        return func(x)
    else:
        return math.factorial(10)/math.factorial(10 - order) * (x **(10 - order))   

def maximum(xArrayGraph, n, freq):
    a = [0]*freq
    for i in range(freq):
        a[i] = abs(derivative(xArrayGraph[i], n))
    return max(a)

n = 20
left = -1
right = 1.5
xArray = []
yArray = []

print("Таблица значений")
j = 0
print("x:          y:")
for i in np.linspace(left, right, n):
    xArray.append(i)
    #xArray.append(Chebushev(left, right, j, n))
    yArray.append(func(xArray[j]))
    print("%-5f   %7f" % (xArray[j],yArray[j]))
    j += 1

freq = 100
xArrayGraph = [0]*freq
yArrayGraph = [0]*freq

print("Полином Лагранжа: ")
yArrayLag = [0]*freq
for i in range(freq):
    xArrayGraph[i] = left + (right-left)*i/(freq-1)
    yArrayGraph[i] = func(xArrayGraph[i])
    yArrayLag[i] = lagrange(xArrayGraph[i], n, xArray, yArray)

#plt.plot(xArrayGraph, yArrayGraph)
#plt.plot(xArrayGraph, yArrayLag)

print("Полином Ньютона: ")
yArrayNew = [0]*freq
for i in range(freq):
    yArrayNew[i] = newton(xArrayGraph[i], n, xArray, yArray)
    
#plt.plot(xArrayGraph, yArrayGraph)
#plt.plot(xArrayGraph, yArrayNew)

print("Погрешности: ")
yArraySin = [0]*freq
sinLag = [0]*freq
sinNew = [0]*freq
for i in range(freq):
    yArraySin[i] = sin(xArrayGraph[i], n, xArray, xArrayGraph, freq)
    sinLag[i] = abs(yArrayGraph[i] - yArrayLag[i])
    sinNew[i] = abs(yArrayGraph[i] - yArrayNew[i])

plt.plot(xArrayGraph, yArrayGraph)
plt.plot(xArrayGraph, yArraySin)
plt.plot(xArrayGraph, sinLag)
plt.plot(xArrayGraph, sinNew)

plt.show()

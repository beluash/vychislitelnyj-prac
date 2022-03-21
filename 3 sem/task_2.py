import math
import numpy as np
import matplotlib.pyplot as plt

#интервал
a = -2
b = 2
n = 8
h = (b-a)/n

#точное значение интеграла
integral = 2.83229367309428

def func(x):
    return abs(math.sin(x))

#Считаем значения интеграла
def midrectangles(h, n):
    left = a
    result = 0
    for i in range(n):
        result += h * func(left + h/2)
        left += h
    return result

def trapezes(h, n):
    left = a
    result = 0
    for i in range(n):
        result += h/2 * (func(left) + func(left + h))
        left += h
    return result

def simpson(h, n):
    left = a
    result = 0
    for i in range(n):
        result += h/6 * (func(left) + 4*func(left + h/2) + func(left + h))
        left += h
    return result

def gauss2(h, n):
    left = a
    result = 0
    for i in range(n):
        x1 = left + h/2 + (3**0.5)/6 * h
        x2 = left + h/2 - (3**0.5)/6 * h
        result += h/2 * (func(x1) + func(x2))
        left += h
    return result

def gauss3(h, n):
    left = a
    result = 0
    for i in range(n):
        x1 = left + h/2 + h/2 * ((3/5)**0.5)
        x2 = left + h/2
        x3 = left + h/2 - h/2 * ((3/5)**0.5)
        result += h/2 * ((5/9)*func(x1) + (8/9)*func(x2) + (5/9)*func(x3))
        left += h
    return result

#Считаем практический порядок точности
def torder(f1, f2):
    e1 = integral - f1
    e2 = integral - f2
    return (math.log(e1/e2, 2))

#Считаем Рунге
def runge(f1, f2, ord):
    return abs(f2-f1)/(2**ord - 1)

xGraph = ["Средние прямоугольники", "Трапеции", "Симпсон", "Гаусс с двумя узлами", "Гаусс с тремя узлами"]
yGraph = [abs(integral - midrectangles(h, n)), abs(integral - trapezes(h, n)), abs(integral - simpson(h, n)), abs(integral - gauss2(h, n)), abs(integral - gauss3(h, n))]

print("f(x) = |sin(x)|")
print("n =", n)
print("Точное значение интеграла:", integral)
print("Формула средних прямоугольников:", midrectangles(h, n))
print("Формула трапеций:", trapezes(h, n))
print("Формула Симпсона:", simpson(h, n))
print("Формула Гаусса с двумя узлами:", gauss2(h, n))
print("Формула Гаусса с тремя узлами:", gauss3(h, n))
print("Экспериментальный порядок точности для средних прямоугольников:", torder(midrectangles(h, n), midrectangles(h/2, n*2)))
print("Экспериментальный порядок точности для трапеций:", torder(trapezes(h, n), trapezes(h/2, n*2)))
print("Экспериментальный порядок точности для Симпсона:", torder(simpson(h, n), simpson(h/2, n*2)))
print("Экспериментальный порядок точности для Гаусса 2:", torder(gauss2(h, n), gauss2(h/2, n*2)))
print("Экспериментальный порядок точности для Гаусса 3:", torder(gauss3(h, n), gauss3(h/2, n*2)))
print("Рунге для прямоугольников:", runge(midrectangles(h, n), midrectangles(h/2, n*2), 2))
print("Рунге для трапеций:", runge(trapezes(h, n), trapezes(h/2, n*2), 2))
print("Рунге для Симпсона:", runge(simpson(h, n), simpson(h/2, n*2), 4))
print("Рунге для Гаусса 2:", runge(gauss2(h, n), gauss2(h/2, n*2), 4))
print("Рунге для Гаусса 3:", runge(gauss3(h, n), gauss3(h/2, n*2), 6))

#График погрешностей
plt.plot(xGraph, yGraph)
plt.ylabel('Значение погрешности')
plt.title("Погрешности")
plt.show()

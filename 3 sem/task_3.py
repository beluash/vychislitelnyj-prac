import math
import numpy as np
import matplotlib.pyplot as plt

left1 = -3
right1 = -2.5
left2 = -0.25
left2_ = -0.3
right2 = 0.25
left3 = 2.5
right3 = 3
eps = 10 ** (-6)
# teta = 2/(m1 + M1) q = (M1 - m1)/(m1 + M1)
teta = 0.35184
q = 0.071756

def func(x):
    return x**2+4*math.cos(x)-4

def dfunc(x):
    return 2*x-4*math.sin(x)

def ddfunc(x):
    return 2-4*math.cos(x)

def method1(a, b):
    count = 0
    while (abs(b-a) > eps):
        count += 1
        c = (a+b)/2
        if func(c) == 0: return c, count
        if (func(a)*func(c) < 0):
            b = c
        else: 
            a = c
    return c, count

def method2(a, b):
    count = 0
    x1 = b
    x0 = a
    while (eps < abs(func(x0)/13.936)):
        x0 = x1
        count += 1
        if func(x1) == 0: return x1, count
        x1 = x0 - func(x0)/dfunc(x0)
        if ((x1 <= a) or (x1 >= b)): 
            x1 = (x0+x1)/2
    return x1, count

def method3 (a, b):
    count = 0
    x0 = b
    x1 = a
    while (abs(x1-x0) > eps):
        count += 1
        if func(x1) == 0: return x1, count
        x0 = x1 - (x1 - x0)*func(x1)/(func(x1) - func(x0))
        x1 = x0 - (x0 - x1)*func(x0)/(func(x0) - func(x1))
    return x1, count

def method4(a, b):
    count = 0
    x0 = a
    x1 = b
    temp = x0
    y = temp
    flag = 1
    while (flag):
        count += 1
        y = temp
        temp = x1 - (func(x1)*(x1-x0))/(func(x1)-func(x0))
        x0 = x1
        x1 = temp
        if (abs(x1-x0) < eps): flag = 0
    return x1, count

def method5 (a, b):
    count = 0
    x0 = a
    x1 = b
    while (abs(x1-x0) > abs(func(x0)/13.936)):
        count += 1
        if func(x1) == 0: return x1, count
        x0 = x1
        x1 = x0 - func(x0)/dfunc(x0) - (ddfunc(x0)*(func(x0)**2))/(2*(dfunc(x0)**3))
    return x1, count

def method6(a, b):
    count = 0
    x0, x2 = a, b
    e, result = 1, 0
    while (abs(e-result) > eps):
        count += 1
        e = result
        x1 = (x0 + x2) / 2
        A = ((func(x2) - func(x1))/(x2 - x1) - (func(x1)-func(x0))/(x1 - x0)) / (x2 - x0)
        B = (func(x2) - func(x1))/(x2 - x1) + (x2 - x1)*A
        C = func(x2)
        root1 = x2 + (-B + (B*B - 4*A*C)**0.5) / 2 / A
        root2 = x2 + (-B - (B*B - 4*A*C)**0.5) / 2 / A
        if ((root1 >= x0 and root1 <= x2) or (root1 >= x2 and root1 <= x0)): result = root1
        else: result = root2
        if (func(x0)*func(result) <  0): x2 = result
        else: x0 = result
    return result, count

def method7 (a, b):
    count = 0
    x0, x1 = a, b
    n = (int)((math.log(eps/abs(b-a)) // math.log(q)) + 1)
    for i in range(n):
        count += 1
        x1 = x0 + teta*func(x0)
        x0 = x1
    return x1, count


print("roots: 0, 2.7831147565030203006, -2.7831147565030203006")
print("Метод дихотомии: ")
print(method1(left1, right1), method1(left2, right2), method1(left3, right3))
print("Метод Ньютона: ")
print(method2(left1, right1), method2(left2, right2), method2(left3, right3))
print("Метод хорд: ")
print(method3(left1, right1), method3(left2_, right2), method3(left3, right3))
print("Метод секущих: ")
print(method4(left1, right1), method4(left2_, right2), method4(left3, right3))
print("Метод Чебышева: ")
print(method5(left1, right1), method5(left2, right2), method5(left3, right3))
print("Метод парабол: ")
print(method6(left1, right1), method6(left2, right2), method6(left3, right3))
print("Метод простой итерации: ")
print(method7(left1, right1), method7(left2, right2), method7(left3, right3))

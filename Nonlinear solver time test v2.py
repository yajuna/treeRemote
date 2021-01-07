import numpy as np
import scipy.optimize
import time
from functools import partial


def F(T):
    sigma = 5.67*10**-8
    T2 = 298
    return sigma*(T2**4-T**4)+30


def linear(t):
    m = 1
    b = 3
    return m*t+b


def quadratic(t):
    a = -1 # making sure it has two solutions
    b = 1
    c = 1
    return a*t**2+b*t+c


sol = (30/5.67*10**(-8)+298**4)**(1/4)
print("sol = ")
print(sol)

# IG = Initial Guess
IG = 10**6

# Method 1 - optimize.root
tic = time.perf_counter()
sol1 = scipy.optimize.root(F, IG)
toc = time.perf_counter()
time_method1 = toc - tic
print("sol1 = ")
print(sol1)
print('time method 1 = ')
print(time_method1)

# Method 2 - optimize.bryoden1
tic = time.perf_counter()
sol2 = scipy.optimize.broyden1(F, IG)
toc = time.perf_counter()
time_method2 = toc - tic
print("sol2 = ")
print(sol2)
print("time method 2 = ")
print(time_method2)


# Method 2 shown to be faster, so testing that on additional functions

# Linear test using bryoden1 (solution should be -3)
IG = 0
tic = time.perf_counter()
sol3 = scipy.optimize.broyden1(linear, IG)
toc = time.perf_counter()
time_method3 = toc - tic
print("sol3 = ")
print(sol3)
print("time method linear = ")
print(time_method3)

# Quadratic test using bryoden1 (solutions are -0.618 and 1.618)
# Shows that bryoden1 only finds nearest solution
tic = time.perf_counter()
sol4 = scipy.optimize.broyden1(quadratic, IG)
toc = time.perf_counter()
time_method4 = toc - tic
print("sol4 = ")
print(sol4)
print("time method quadratic = ")
print(time_method4)

# Quadratic test 2 using bryoden1
# For some reason this method consistantly finds solution about 10x faster than other quadratic solution
# showing that time varies a lot using this method
IG = 1
tic = time.perf_counter()
sol5 = scipy.optimize.broyden1(quadratic, IG)
toc = time.perf_counter()
time_method5 = toc - tic
print("sol5 = ")
print(sol5)
print("time method quadratic = ")
print(time_method5)

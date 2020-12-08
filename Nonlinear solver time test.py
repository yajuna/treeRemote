import numpy as np
import scipy.optimize
import time
from functools import partial


def F(T):
    sigma = 5.67*10**-8
    T2 = 298
    return sigma*(T2**4-T**4)+30


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

# Method 3 - Selina's code ****needs work***
def findtemp_on_array(T_array, u_array):
    print('T is ' + str(T_array))
    print('u is ' + str(u_array))
    ans = [F2(T, u) for T, u in zip(T_array, u_array)]
    print('result array is ' + str(ans))
    return ans  # this line of code calls findtemp on each corresponding pair of u, T values (tuple)


def F2(T, u):
    sigma = 1.38*10**-23
    return sigma*(u**4-T**4)+30


#sol3 = findtemp_on_array(IG, 2*IG)
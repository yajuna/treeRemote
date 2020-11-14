import numpy as np
import scipy.optimize
from functools import partial


def F(T):
    sigma = 1.38*10**-23
    T2 = 298
    return sigma*(T2**4-T**4)-30


sol = scipy.optimize.root(F, 10**7)
print(sol)
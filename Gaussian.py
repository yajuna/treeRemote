# 6/2/2020
# Test source terms to help debug the code
# Michael Hockman

import numpy as np

# Setting up the variables, assuming we are starting at midnight the previous night
variance = 1
mean = 12  # peak temp at noon
t_end = 24  # Assuming we start at midnight and end at midnight
size = 100  # How many source terms we are going to get
t = np.linspace(0, 24, size)
A = 75


# takes in a time ti, variance v, mean m, and amplitude A
# returns the output of a Gaussian function with mean m and variance v
def sourceTerms(ti, A, m, v):
    coefficient = A / (np.sqrt(2 * np.pi) * v)
    exponential = -(ti - m) ** 2 / (2 * v ** 2)
    return coefficient * np.exp(exponential)


terms = []  # Source terms
segments = []  # Time intervals each source term is mapped to, in 24 hour clock time
for i in range(0, size):
    segments.append([[24 / size * i], [24 / size * (i + 1)]])
    terms.append(sourceTerms(i, A, mean, variance))


print(segments)
print(terms)

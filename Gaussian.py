# 6/2/2020
# Test source terms to help debug the code
# Michael Hockman

import numpy as np
import matplotlib.pyplot as plt

# Setting up the variables, assuming we are starting at midnight the previous night
variance = 2
mean = 12  # peak temp at noon
t_end = 24  # Assuming we start at midnight and end at midnight
size = 24  # How many source terms we are going to get, for testing purposes
t = np.linspace(0, 24, size)  # time starting at hour 0 and ending at hour 24
A = 75  # variable for the Gaussian, controls max value
B = 0  # variable for the Gaussian, controls the offset

# takes in a time ti
# returns the output of a Gaussian function with mean and variance depending on the variables defined above
def sourceTerms(ti):
    coefficient = A / (np.sqrt(2 * np.pi) * variance)
    exponential = -(ti - mean) ** 2 / (2 * variance ** 2)
    return B + coefficient * np.exp(exponential)


terms = []  # Source terms
segments = []  # Time intervals each source term is mapped to, in 24 hour clock time
for i in range(0, size):
    segments.append([[24 / size * i], [24 / size * (i + 1)]])
    terms.append(sourceTerms(i))

# print(segments)
# print(terms)

plt.xlabel("time (hrs)")
plt.ylabel("temperature (c)")
plt.plot(t, terms)
plt.show()
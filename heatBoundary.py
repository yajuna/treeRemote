# 4/13/2020
# Boundary conditions and shadow tracking for tree
# Michael Hockman

import numpy as np

# size = how many subsections each half of the tree is split into
size = 50
w = np.pi / (size - 1)
theta0 = -np.pi / 2  # setting it so it starts south, so sun is shining on the east side
s = []


# R = radius of tree, omega = speed the shaded/non-shaded halves are rotating, t0 = time of day
# returns a set of [x, y] points on the boundary of tree
def pointP(R, omega, t0):
    ## Clears P ev
    p = []
    x0 = R * np.cos(theta0 + omega * t0)
    if abs(x0) < 10 ** -10:
        x0 = 0
    y0 = R * np.sin(theta0 + omega * t0)
    if abs(y0) < 10 ** -10:
        y0 = 0
    # first point
    p.append([x0, y0])
    x1 = R * np.cos(theta0 + omega * t0 + np.pi)
    if abs(x1) < 10 ** -10:
        x1 = 0
    y1 = R * np.sin(theta0 + omega * t0 + np.pi)
    if abs(y1) < 10 ** -10:
        y1 = 0
    # second point
    p.append([x1, y1])
    return p


# creates a set of sets, [[[x(theta0), y(theta0)], [x(theta0+pi), y(theta0+pi)]],  [x(theta1), y(theta1)],
# [x(theta1+pi), y(theta1+pi)]] ....] continues up to theta_(size-1)
for i in range(size):
    s.append(pointP(1, w, i))

# s = sum from i = 0 to i = size - 1 [[r*cos(theta+w*i), r*sin(theta+w*i)], [r*cos(theta+w*i+pi), r*sin(theta+w*i+pi)]]
# can remove this when we're done testing this
print(s)
# split time into a 12 hour period so it ends at midnight, assuming it begins at noon
t = np.linspace(0, 12, size)


#  https://www.desmos.com/calculator/4kh1zxm6fl <--- graph with adjustable parameters 
# omega0 controls quickly the arctan portion goes to pi/2
# omega1 controls how much the trees temperature varies in time over the day
# A is the nighttime temperature
# B controls how much hotter than sunnyside is than the shaded side
# C controls how fast the tree approaches its nighttime temperature where Tsun = Tshade
# C should be approx .3 < C < 1 for the temp to converge around
# ti is the time at a certain instant
def sunnyside(A, B, omega0, omega1, ti, C):
    return (A + B * np.exp(-C * ti)) - (1 + np.exp(-C * ti) * np.cos(omega1 * ti)) * np.arctan(omega0 * ti)


# omega0 controls quickly the arctan portion goes to pi/2
# omega1 controls how much the trees temperature varies in time over the day
# A is the nighttime temperature
# C controls how fast the tree approaches its nighttime temperature where Tsun = Tshade
# approx 0 < C < 1
# ti is the time at a certain instant
# n controls how faster the shaded side reaches its asymptotic minimum, essentially controls how much colder
# the shaded side is than the sunny side,
# 1 < n, or else tshade > tsun sometimes
def shadyside(A, n, omega0, omega1, ti, C):
    return A - (1 + np.exp(-C * ti) * np.cos(omega1 * ti)) * np.arctan(n * omega0)


# temp sun, temp shade
Tsun = []
Tshade = []
Ttotal = []
# fills up vectors Tsun and Tshade
#     [0, 1, 2, 3,  4,      5]
# X = [A, B, C, n, omega0, omega1]
X = [4, 6, 0.48, 1.4, 0.6, 1.5]

# sunnyside(A, B, omega0, omega1, ti, C)
# shadyside(A, n, omega0, omega1, ti, C)
# assigning a temp to each side a constant temperature value for every time defined on line 47
for i in range(0, size):
    Tsun.append(sunnyside(X[0], X[1], X[4], X[5], t[i], X[2]))
    Tshade.append(shadyside(X[0], X[3], X[4], X[5], t[i], X[2]))

# Bringing the two halves together
Ttotal = Tsun + Tshade

# printing the output for testing purposes
print("Tshade = " + str(Tshade))
print("Tsun = " + str(Tsun))
# Ttotal[0,size-1] = Tsun
# Ttotal[size, 2size-1] = Tshade
print("Ttotal = " + str(Ttotal))

# should = 0, since there are 2 coordinates in "s" for every 1 component there is in "Ttotal" (temp in that interval)
print(np.size(Ttotal) - 1 / 2 * np.size(s))

# look at main code or heat equation paper (potter & anderson) equation 10 to create the internal temp

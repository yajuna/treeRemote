# 4/13/2020
# Boundary conditions and shadow tracking for tree
# Michael Hockman

import numpy as np

# size = how many subsections each half of the tree is split into
size = 10
w = np.pi / (size - 1)
theta0 = -np.pi / 2  # setting it so it starts south, so sun is shining on the east side
s = []


# R = radius of tree, omega = speed the shaded/non-shaded halves are rotating, t0 = time of day
# returns a set of [x, y] points on the boundary of tree
def pointP(R, omega, t0):
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

# prints sum from i = 0 to i = size-1 [[r*cos(theta)+w*i, r*sin(theta+w*i)], [r*cos(theta+w*i+pi), r*sin(theta+w*i+pi)]]
# can remove this when we're done testing this
print(s)

# since average of a sinusoid from 0 to 2pi is 0, we will use average from 0 to pi/2
N = size
t = np.linspace(0, np.pi / 4, N)


# depending on what units of heat we're using changing 30 to something more reasonable
def sunnyside(A1, omega, t):
    return np.average(30 + A1 * np.cos(omega * t))


def shadyside(A0, omega, t):
    return 0.8 * sunnyside(A0, omega, t)


# temp sun, temp shade
Tsun = []
Tshade = []

# fills up vectors Tsun and Tshade
for i in range(0, N):
    Tsun.append(sunnyside(1, 1, t[i]))
    Tshade.append(0.8 * Tsun[i])


# assume temp inside the tree is constant with theta for now
# omega0 = omega used for sunnyside and shadyside functions
# omega1 = the angular velocity of the temperature inside the tree with respect to radius
# Makes a sinusoid with magnitude 1 / 2 * (shadyside(A, w0, time) + sunnyside(A, w0, time))
def internal_temp(A, w0, w1, r, time):
    return 1 / 2 * (shadyside(A, w0, time) + sunnyside(A, w0, time)) * cos(w1 * r)

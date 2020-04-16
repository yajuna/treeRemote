#4/13/2020
#Boundary conditions and shadow tracking for tree
#Michael Hockman

import numpy as np

size = 10
w = np.pi / (size-1)
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

# creates a set of sets, [[[x(theta0), y(theta0)], [x(theta0+pi), y(theta0+pi)]],  [x(theta1), y(theta1)], [x(theta1+pi), y(theta1+pi)]] ....]
# continues up to theta_(size-1)
for i in range(size):
    s.append(pointP(1, w, i))

print(s)

t = np.linspace(0,5,500)
def sunnySide(A1, omega, t):
    return np.average(A1*np.cos(omega*t))



def shadySide(A0, omega, t):
    return 0.8*sunnySide(A0, omega, 0, t)
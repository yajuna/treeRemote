#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:27:41 2020

modify heat2D to polar form. Add Gausian source term
ref: heat2D.py

@author: yajun
"""

import numpy as np
import matplotlib.pyplot as plt

# tree trunk size, mm
angle = 10
r = 10.
# intervals in r-, phi- directions, mm
dr = dphi = 0.1
# Thermal diffusivity of steel, mm2.s-1
D = 4.

Tcool, Thot = 300, 700

nr, nphi = int(r/dr), int(angle/dphi)
dr2, dphi2 = dr*dr, dphi*dphi
dt = dr2 * dphi2 / (2 * D * (dr2 + dphi2))

u0 = Tcool * np.ones((nr, nphi))
u = u0.copy()

# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)
r, cx, cy = 2, 5, 5
r2 = r**2
for i in range(nr):
    for j in range(nphi):
        p2 = (i*dr-cx)**2 + (j*dphi-cy)**2
        if p2 < r2:
            u0[i,j] = Thot

def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dr2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dphi2 )

    u0 = u.copy()
    return u0, u


# Number of timesteps
nsteps = 1001#101
# Output 4 figures at these timesteps
mfig = [0, 10, 50, 100]#[700, 800, 900, 1000]
fignum = 0
fig = plt.figure()
for m in range(nsteps):
    u0, u = do_timestep(u0, u)  + np.exp((0.313-dt* m)**2) # modify to add Gaussian source
    print(np.max(u))
    if m in mfig:
        fignum += 1
        print(m, fignum, np.max(u))
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tcool,vmax=Thot)
        ax.set_axis_off()
        ax.set_title('{:.1f} ms'.format(m*dt*1000))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('$T$ / K', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.show()
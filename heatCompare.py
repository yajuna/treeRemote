#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 13:42:04 2020

solve the heat equation in polar with FD. see c14_2d_disk_heat.pdf for ref
ut = k(urr + ur/r + uthetatheta/r^2), 0<r<1, 0<theta<2pi
u(1,theta,t) = 0, 0<theta<2pi, outer bdry condition
u(r,0,t) = u(r,2pi,t),  0<r<1
u(r,theta,t) bounded
u(r,theta,0) = f(r,theta) = (r-r^3)sin(theta)

@author: yajun
"""


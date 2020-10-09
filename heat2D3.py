#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 13:21:57 2020

@author: yajun
"""

# In[1]:

get_ipython().magic(u'pylab inline')


# In[3]:

import scipy as sp
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
 

class Heat_Equation(object):
    """
    Class which implements a numerical solution of the 2d heat equation
    """
    def __init__(self, dx, dy, a, kind, timesteps = 1):
                 self.dx = dx # Interval size in x-direction.
                 self.dy = dy # Interval size in y-direction.
                 self.a = a # Diffusion constant.
                 self.timesteps = timesteps  #Number of time-steps to evolve system.
                 self.dx2 = dx**2
                 self.dy2 = dy**2
                 self.nx = int(1/dx)
                 self.ny = int(1/dy)
                # For stability, this is the largest interval possible
                # for the size of the time-step:
                 self.dt = self.dx2*self.dy2/( 2*a*(self.dx2+self.dy2) )
                 self.u,self.ui = self.get_initial_conditions(kind)
                    
    def get_initial_conditions(self, kind):
        # Start u and ui off as zero matrices:
        ui = sp.zeros([self.nx,self.ny])
        u = sp.zeros([self.nx,self.ny])
        # Now, set the initial conditions (ui).
        for i in range(self.nx):
            for j in range(self.ny):
                if kind == "two_circles":
                    p = (i*self.dx - 0.5)**2 + (j*self.dy - 0.5)**2
                    if ( p  <= .03 and p >= 0.020 ):
                        ui[i,j] = 1
                    elif ( p  <= .108 and p >= 0.09 ):
                        ui[i,j] = 1
                elif kind == "part_circle":
                    p = (i*self.dx - 0.5)**2 + (j*self.dy - 0.5)**2
                    if ( p  <= .01 and p >= 0.009):
                        ui[i,j] = 1
                elif kind == "four_blobs":
                    p = (i*self.dx - .4)**2 + (j*self.dy - .4)**2
                    if ( p  <= 0.02 and p >= 0.02):
                        ui[i,j] = 1
                elif kind == "two_blobs":
                    p = (i*self.dx - .4)**2 + (j*self.dy - .4)**2
                    if ( p  <= 0.05 and p >= .05):
                        ui[i,j] = 1
                elif kind == "half_moon":
                    p = (i*self.dx-.2)**2+(j*self.dy-.2)**2
                    if ( p <= 0.1 and p >=.05 ):
                        ui[i,j] = 1
                elif kind == "2_lines":
                    p = (i*self.dx-.5)**2+(j*self.dy-j*.4)**2
                    if ( p <= 0.1 and p >=.05 ):
                        ui[i,j] = 1
                elif kind == "circle":
                    p = (i*self.dx-0.5)**2+(j*self.dy-0.5)**2
                    if ( p <= 0.1 and p >=.05 ):
                        ui[i,j] = 1
                elif kind == "lower_blob":
                    p = (i*self.dx-0.)**2+(j*self.dy-0.)**2 
                    if ( p <= 0.1 and p >=.05 ):
                        ui[i,j] = 1
                elif kind  == "two_half_moons":
                    p = (i*self.dx - 0.5)**3 + (j*self.dy - 0.5)**3
                    if ( p  <= .03 and p >= 0.020 ):
                        ui[i,j] = 1
                    elif ( p  <= .108 and p >= 0.09 ):
                        ui[i,j] = 1
        return u,ui
    
    def evolve_ts(self):
        self.u[1:-1, 1:-1] = self.ui[1:-1, 1:-1] + self.a*self.dt*( (self.ui[2:, 1:-1] - 2*self.ui[1:-1, 1:-1] + self.ui[:-2, 1:-1])/self.dx2 + (self.ui[1:-1, 2:] - 2*self.ui[1:-1, 1:-1] + self.ui[1:-1, :-2])/self.dy2 )
        self.ui = self.u.copy()
        
def evolve_ts(u, ui):
    global nx, ny
    """
    This function uses two plain Python loops to
    evaluate the derivatives in the Laplacian, and
    calculates u[i,j] based on ui[i,j].
    """
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            uxx = ( ui[i+1,j] - 2*ui[i,j] + ui[i-1, j] )/ dx2
            uyy = ( ui[i,j+1] - 2*ui[i,j] + ui[i, j-1] )/ dy2
            u[i,j] = ui[i,j]+dt*a*(uxx+uyy)


# In[4]:

from tempfile import NamedTemporaryFile

def anim_to_mp4(anim,k):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            newname = "/tmp/" + k + f.name.split("/")[-1]
            print(newname)

            anim.save(newname, fps=20, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])    
    return None

def save_animation(anim,k):
    plt.close(anim._fig)
    return anim_to_mp4(anim,k)


# In[5]:

from matplotlib import animation
#initial_shapes = ["two_circles", "part_circle", "four_blobs", "two_blobs", "half_moon", "two_half_moons", "2_lines", "circle", "lower_blob"]
initial_shapes = ["two_circles"]

for k in initial_shapes:
    test_heat = Heat_Equation(0.001,0.001,.5,k,1)    
    
    # First set up the figure, the axis, and the plot element we want to animate
    
    
    fig = plt.figure()
    img = plt.subplot(111)
    
    im = img.imshow(test_heat.ui, cmap=cm.get_cmap("hot"), interpolation='nearest', origin='lower')
    im.figure = fig
    fig.colorbar(im)
                
    def animate(i,im):
        if i % 50 == 0:
            print(i)
        test_heat.evolve_ts()
        im.set_array(test_heat.ui)

        return [im]
        
        
    anim = animation.FuncAnimation(fig, animate, frames=2000, fargs=(im,),  interval=30, blit=True)
    save_animation(anim,k)

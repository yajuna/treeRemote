#  Created 7/14/2020 by Michael Hockman
#  First attempt at combining the source terms into one term
#  To do this I'm using figure 4 from the potter anderson paper to find the conduction "k",
#  the convective heat exchange "H", the insolation "Sdir+Sdif", and net infrared "IRin - IRout"

import matplotlib.pyplot as plt
import numpy as np

albedo = .3  # "For albedo, we assumed a value of 0.3" - Potter Anderson page 3

#  North variables
insolationN = [28, 30, 29, 27, 26, 25, 26, 27, 28, 30, 28, 14, 0]  # IR_in-IR_out
netInfN = [13, -9, -7, -2, -5, -7, -8, -6, -7, -14, -16, -28, -8]  # S_dir+S_dif
convectionN = [13, -9, -7, -2, -5, -7, -8, -9, -16, -18, -33, -30, -22]  # H
# conductionN = [15, -10, -8, -11, -8, -3, -1, -2, 1, 7, 20, 22, 20]  # k
totalN = []

#  Make sure t = [0, 1, 2, ... , 12]
#  because we need 13 on the t-axis for the 13 data points the Potter Anderson paper provides.
t = np.linspace(0, 12, 13)

# x is the number of points the linear interpolation will have
# currently set to be 10x the size of t
x = np.linspace(0, 12, 1000)

#  interpolate the north data
#  Also turn into an array so it sums how I want it to for totalN
insolationNI = np.array(np.interp(x, t, insolationN))
netInfNI = (1-albedo)*np.array(np.interp(x, t, netInfN))  # Scaled by (1-alpha)
convectionNI = np.array(np.interp(x, t, convectionN))

# Sum all the lines up
totalN = insolationNI + netInfNI + convectionNI


#  South variables
insolationS = [20, 180, 240, 330, 370, 390, 370, 330, 240, 180, 40, 30, 20]  # IR_in-IR_out
netInfS = [0, -10, -15, -20, -25, -30, -25, -20, -18, -16, -14, -12, -6]  # S_dir+S_dif
convectionS = [0, -10, -25, -30, -35, -40, -35, -30, -28, -22, -18, -12, -6]  # H
# conductionS = [0, -100, -150, -190, -195, -200, -180, -130, -80, -30, -28, -10, -6]  # k
totalS = []

#  interpolate the south data
#  Also turn into an array so it sums how I want it to for totalS
insolationSI = np.array(np.interp(x, t, insolationS))
netInfSI = (1-albedo)*np.array(np.interp(x, t, netInfS))  # Scaled by (1-alpha)
convectionSI = np.array(np.interp(x, t, convectionS))

# Sum all the lines up
totalS = insolationSI + netInfSI + convectionSI



print(totalS)
print(type(totalS))
print(np.size(totalS))



#  Uncomment if you want to see the plot
plt.plot(x[:], totalS[:], "-g", label="south")
plt.plot(x[:], totalN[:], "-b", label="north")
plt.axis([0,12,-60,350])
plt.legend()
plt.xlabel("Time since 7:00am (hrs)")
plt.ylabel("Energy flux (W/m^2)")
plt.grid()
plt.show()


#  input: t0, the time at which the source term will be evaluated at
#  return: North, the total source terms at time t0
#  if it can't be evaluated at t0 exactly, it will find the closest available point to evaluate
def sourceTermsN(t0):
    t0 = int(t0 / 0.008)
    return totalN[t0]



#  input: t0, the time at which the source term will be evaluated at
#  return: S evaluated at t0, the total source terms at time t0
#  if it can't be evaluated at t0 exactly, it will find the closest available point to evaluate
def sourceTermsS(t0):
    t0 = int(t0/0.008)
    print(t0)
    return totalS[t0]


# Simply returns an array of size t instead of a single element at point t0
# t is an input to make it match with the Heat1D file
def sourceTermsArrayS(t):
    return totalS[np.arange(1, np.size(t))]


# Simply returns an array of size t instead of a single element at point t0
# t is an input to make it match with the Heat1D file
def sourceTermsArrayN(t):
    return totalN[np.arange(1, np.size(t))]

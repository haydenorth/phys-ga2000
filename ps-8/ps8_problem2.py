# Graduate Computational Physics
# Problem Set 8
# Problem 2
# By Hayden Orth

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

# constants
global sigma, r, b
sigma = 10
r = 28
b = 8/3

# function for solving 3 diff eq
def f(t, r0):
    x, y, z = r0
    fx = sigma * (y - x)
    fy = r*x - y - x*z
    fz = x*y - b*z
    return np.array([fx, fy, fz])

# function to plot y vs t
def ploty_t(y, t):
    plt.plot(t, y)
    plt.title('y vs. time')
    plt.xlabel('Time (s)')
    plt.ylabel('y')
    plt.savefig('ps8_p2_yt.png')

# plot z vs x
def plotz_x(z, x):
    plt.plot(x, z)
    plt.title('Strange Attractor')
    plt.xlabel('x')
    plt.ylabel('z')
    plt.savefig('ps8_p2_zx.png')

# solve diff eqs with scipy
t = np.linspace(0,50,1200)
sol = integrate.solve_ivp(f, (0, 50), (0, 1, 0), method ='Radau',
                          t_eval = t, rtol=1e-6)
x, y, z = sol.y

# make plots
#ploty_t(y, t)
plotz_x(z, x)
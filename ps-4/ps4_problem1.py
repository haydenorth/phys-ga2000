# Graduate Computational Physics
# Problem Set 4
# Problem 1
# By Hayden Orth

from gaussxw import gaussxwab
import numpy as np
import matplotlib.pyplot as plt

# integrand of Cv integral
def integrand(x):
    return (x**4 * np.exp(x)) / (np.exp(x) - 1)**2

# function for calculating integral via Gaussian quadrature
def gquad(a, b, N):
    x, w = gaussxwab(N,a,b)
    s = sum(integrand(x)*w)
    return s

# function for calculating Cv for given T, N
def cv(T, N):
    upper_bound = 428 / T
    integral = gquad(0, upper_bound, N)
    constants = 9*0.001*6.022e28*1.38e-23* (T / 428)**3
    return constants*integral

N = 50

# gather Cv data for T from 5 K to 500 K
Tpoints = np.arange(5, 505, 5)
cvpoints = []
for T in Tpoints:
    cvpoints.append(cv(T, N))

# create plot for part b
plt.plot(Tpoints, cvpoints)
plt.title('Heat Capacity of a Solid')
plt.xlabel('Temperature (K)')
plt.ylabel('Heat Capacity (J/K)')
plt.savefig('heatcapacity.png')

# test convergence for T = 500 K
T = 500
Npoints = np.arange(10, 80, 10)
cvpoints2 = []
for N in Npoints:
    cvpoints2.append(cv(T, N))

# create plot for part c
plt.plot(Npoints, cvpoints2)
plt.title('Cv Calculation Convergence')
plt.xlabel('Number of Sample Points')
plt.ylabel('Heat Capacity (J/K)')
plt.savefig('cv_convergence.png')

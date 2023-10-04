# Graduate Computational Physics
# Problem Set 4
# Problem 3
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt
import math
from gaussxw import gaussxwab

# recursive function for Hermite polynomials, takes integer n and an array x, returns Hn(x)
def H(n, x):
    if n==0:
        return np.ones(x.shape)
    elif n==1:
        return 2*x
    else:
        return 2*x*H(n-1, x)-2*(n-1)*H(n-2, x)

# function that returns the harmonic oscillator wavefunction for given n and x
def psi(n, x):
    return (1 / np.sqrt(2**n * math.factorial(n) * np.sqrt(np.pi))) * np.exp(-x**2 / 2) * H(n, x)

# plot first four wavefunctions
x = np.linspace(-4, 4, 100)
psi0 = psi(0, x)
psi1 = psi(1, x)
psi2 = psi(2, x)
psi3 = psi(3, x)
# make plot
plt.plot(x, psi0, 'k-', label='n = 0')
plt.plot(x, psi1, 'r-', label='n = 1')
plt.plot(x, psi2, 'g-', label='n = 2')
plt.plot(x, psi3, 'b-', label='n = 3')
plt.legend()
plt.title('Harmonic Oscillator Wavefunctions')
plt.xlabel('x (hbar / mw)^(1/2)')
plt.ylabel('Psi(x)')
plt.savefig('ps4_p3_a.png')

'''
# part b: plot function for n = 30 on x = [-10, 10]
xb = np.linspace(-10, 10, 500)
psi30 = psi(30, xb)
plt.plot(xb, psi30)
plt.title('Harmonic Oscillator: n = 30')
plt.xlabel('x')
plt.ylabel('Psi(x)')
plt.savefig('ps4_p3_b.png')
'''

# part c: calculate rms value of x**2
# integrand of rms integral. includes a change of variables to evaluate integral from -inf to inf
def integrand(x):
    return ((1 + x**2) / (1 - x**2)**2) * (x / (1 - x**2))**2 * np.abs(psi(5, (x / (1-x**2))))**2

# function for calculating integral via Gaussian quadrature
def gquad(a, b, N):
    x, w = gaussxwab(N,a,b)
    s = sum(integrand(x)*w)
    return s
# calculate integral and print rms value
integral = gquad(-1, 1, 100)
print(np.sqrt(integral))

# part d: use Gauss-Hermite quadrature
from scipy.special import roots_hermite

# new function for psi to be used in Gauss-Hermite method. Psi2 does not include the exponential term
def psi2(n, x):
    return (1 / np.sqrt(2**n * math.factorial(n) * np.sqrt(np.pi))) * H(n, x)

# function for integrand in part d. does not contain the change of variables previously used
def integrand2(x):
    return x**2 * np.abs(psi2(5, x))**2

# function for integrating with gaussian-hermite quadrature
def gquad_herm(a, b, N):
    x,w = roots_hermite(N)
    s = sum(integrand2(x)*w)
    return s

# calculate integral and print rms value
int_hermite = gquad_herm(-1, 1, 100)
print(np.sqrt(int_hermite))

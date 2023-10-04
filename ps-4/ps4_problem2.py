# Graduate Computational Physics
# Problem Set 4
# Problem 2
# By Hayden Orth

from gaussxw import gaussxwab
import numpy as np
import matplotlib.pyplot as plt

# integrand of T integral
def integrand(x, a):
    return 1 / (np.sqrt(a**4 - x**4))

# function for calculating integral via Gaussian quadrature
def gquad(a, b, N, amp):
    x, w = gaussxwab(N,a,b)
    s = sum(integrand(x, amp)*w)
    return s

# function that calculated period of the oscillator for a given amplitude a
def get_T(a):
    integral = gquad(0, a, 20, a)
    return np.sqrt(8)*integral

# collect data for T as a funciton of a
apoints = np.linspace(0, 2, 50)
periods = []
for a in apoints:
    periods.append(get_T(a))

# create plot of T vs a
plt.plot(apoints, periods)
plt.title('Anharmonic Oscillator: V(x) = x^4')
plt.xlabel('Oscillation Amplitude (m)')
plt.ylabel('Period (s)')
plt.savefig('ps4_p2.png')
# New York University
# PHYS GA 2000
# Problem Set 1
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt

# create function for gaussian
def Gaussian(x, mu, std):
    # arguments:
    # x is an array of x values
    # mu is the mean of the distribution
    # std is the standard deviation of the distribution

    return (1/(std*np.sqrt(2*np.pi))) * np.exp(-((x-mu)**2) / (2 * std**2))

# create arrays for plot data
x = np.linspace(-10, 10, 1000)
y = Gaussian(x, 0, 3)

# create and save plot
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('A Gaussian')
plt.savefig('gaussian.png')

# having trouble running from terminal: module not found error for numpy

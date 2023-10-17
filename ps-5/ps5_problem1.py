# Graduate Computational Physics
# Problem Set 5
# Problem 1: Newman 5.17
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxwab

# part a: graph the value of the integrand
x0 = np.linspace(0, 5, 100)
def integrand0(x, a):
    return x**(a-1) * np.exp(-x)
# generate data for a = 2, 3, 4
integrand2 = integrand0(x0, 2)
integrand3 = integrand0(x0, 3)
integrand4 = integrand0(x0, 4)
# plot curves
plt.plot(x0, integrand2, label='a = 2')
plt.plot(x0, integrand3, label='a = 3')
plt.plot(x0, integrand4, label='a = 4')
plt.title('Gamma Function Integrand')
plt.xlabel('x')
plt.ylabel('Integrand(x)')
plt.legend()
plt.savefig('ps5_p1_a')

# part e
# integrand of gamma integral with eq (5.69) change of variables, c = (a-1)
def integrand(x, a):
    c = (a-1)
    return (c/(1-x)**2)*np.exp((a-1)*np.log((c*x)/(1-x))-((c*x)/(1-x)))

# function for calculating the gamma function for a given a
# uses gaussian quadrature to evaluate the integral
def gamma(a):
    x, w = gaussxwab(100,0,1)
    return sum(integrand(x, a)*w)

print(gamma(1.5))

# part f
# use gamma function to approximate factorial for integer values of a:
# gamma(a) = (a-1)!
values = [3, 6, 10]
for i in values:
    print(gamma(i))
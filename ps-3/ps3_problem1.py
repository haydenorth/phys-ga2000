# Graduate Computational Physics
# Problem Set 3
# Problem 1
# By Hayden Orth

import numpy as np

# function for the function f(x) = x(x-1)
def func(x):
    return x*(x-1)

# function for the derivative of a function
def deriv(x, d):
    delta = 10**d
    return (func(x+delta) - func(x)) / delta

# calculate derivative using different values for delta
for i in [-2, -4, -6, -8, -10, -12]:
    print(deriv(1, i))
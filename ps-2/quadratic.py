# Graduate Computational Physics
# Problem Set 2
# problem 4
# By Hayden Orth

import numpy as np

# function for finding the roots of a quadratic
def quadratic(a, b, c):
    sqrt_term = np.sqrt(b**2 - 4*a*c)
    # calculate roots according to classic quadratic formula if sqrt_term is to greater precision than 2*a
    if len(str(-b + sqrt_term)) >= len(str(2*a)):
        x1 = (-b + sqrt_term) / (2*a)
        x2 = (-b - sqrt_term) / (2*a)
    else: # otherwise use the modified formula
        x1 = (2*c) / (-b - sqrt_term)
        x2 = (2*c) / (-b + sqrt_term)
    
    return (x1, x2)

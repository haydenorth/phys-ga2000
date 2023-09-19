# Graduate Computational Physics
# Problem Set 2
# problem 2
# By Hayden Orth

import numpy as np
import timeit

# functions to call for timing
# test for loop time
def get_time_for():
    mysetup = 'import numpy as np'
    
    mycode = '''
L = 50
M1 = 0
for i in range(-L, L+1):
    for j in range(-L, L+1):
        for k in range(-L, L+1):
            if i == j == k == 0:
                continue
            M1 += ((-1)**(i+j+k)) / (np.sqrt(i**2 + j**2 + k**2))
'''
    print(timeit.timeit(setup=mysetup,
                    stmt=mycode,
                    number=1))
# test no for loop time
def get_time_nofor():
    mysetup = 'import numpy as np'
    
    mycode = '''
L = 50
x0 = np.arange(-L, L+1)
x, y, z = np.meshgrid(x0, x0, x0)

signs_matrix = (-1)**np.abs(x + y + z)

V0 = np.sqrt(x**2 + y**2 + z**2)
V0 = V0*signs_matrix
V0 = V0[V0 != 0]
V0 = 1 / V0

M2 = np.sum(V0)
'''
    print(timeit.timeit(setup=mysetup,
                    stmt=mycode,
                    number=1))

# L: number of atoms in all 3 direction from origin
L = 50

# M1: Madelung constant to be calculated by for loop
M1 = 0

# calculate M via a for loop, skipping when (i,j,k) = (0,0,0)
for i in range(-L, L+1):
    for j in range(-L, L+1):
        for k in range(-L, L+1):
            if i == j == k == 0:
                continue
            M1 += ((-1)**(i+j+k)) / (np.sqrt(i**2 + j**2 + k**2))

print('M via a for loop:')
print(M1)

# calculate M2, the Madelung constant calculated without a for loop
#M2 = 0

# create a 3D grid of coordinate values with meshgrid
x0 = np.arange(-L, L+1)
x, y, z = np.meshgrid(x0, x0, x0)

# make all indices negative where i+j+k is odd
signs_matrix = (-1)**np.abs(x + y + z)

# calculate each term in the series
V0 = np.sqrt(x**2 + y**2 + z**2)
V0 = V0*signs_matrix
V0 = V0[V0 != 0]
V0 = 1 / V0

# sum all the terms together = M2
M2 = np.sum(V0)

print('\nM without for loop:')
print(M2)

print('\ntime it takes using a for loop:')
get_time_for()
print('\ntime it take without using a for loop:')
get_time_nofor()




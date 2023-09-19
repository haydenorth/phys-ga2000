# Graduate Computational Physics
# Problem Set 2
# problem 3
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt

# N x N grid size
N = 2000
# create 2D of c values: c = x + iy
x0 = np.linspace(-2, 2, N)
x, y = np.meshgrid(x0, x0)

# create grid of c values
c = x + (y * 1j)

# perform I iterations over the c grid
I = 100
z = c
for i in range(1, I):
    z = np.square(z) + c

# set any points not in the set ( |z| <=2 ) equal to zero for plotting
z = np.abs(z)
z[z <= 2] = 0

# plot grid of z values
plt.figure()
plt.imshow(z, cmap = 'gray', extent=[-2,2,-2,2])
plt.title('Mandelbrot Set')
plt.savefig('mandlebrot.png')

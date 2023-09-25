# Graduate Computational Physics
# Problem Set 3
# Problem 4
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt

# constants
tau = 3.053*60
mu = np.log(2) / tau

# make an array of 1000 random numbers form uniform distribution
rng = np.random.default_rng()
uniform_dist_nums = rng.random(1000)

# transform the uniformly distributed random numbers to the nonuniform distribution
t_decay = -(1/mu) * np.log(1 - uniform_dist_nums)
t_decay = np.sort(t_decay)

# see how many atoms havent decayed for each time step, create list of number of Tl atoms
tpoints = np.arange(0, 1750)
NTl = []
for t in tpoints:
    NTl.append(t_decay[t_decay>t].size)

# plot the results
plt.plot(tpoints, NTl)
plt.title('Decay of Tl 208')
plt.xlabel('Time (s)')
plt.ylabel('Number of Tl atoms')
plt.savefig('radiodecay_part2')


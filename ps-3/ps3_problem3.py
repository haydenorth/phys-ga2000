# Graduate Computational Physics
# Problem Set 3
# Problem 3
# By Hayden Orth

import numpy as np
import random
import matplotlib.pyplot as plt

# constants, array/list initialization
NBi = 10000
NTl = 0
NPb = 0
NBi_209 = 0

tau_Bi = 60*46
tau_Tl = 60*2.2
tau_Pb = 60*3.3
tmax = 20000
h = 1.0

tpoints = np.arange(0, tmax)
p_Bi = 1 - 2**(-h/tau_Bi)
p_Tl = 1 - 2**(-h/tau_Tl)
p_Pb = 1 - 2**(-h/tau_Pb)
Bipoints = []
Tlpoints = []
Pbpoints = []
Bi209points = []

# main simulation loop
for t in tpoints:
    Bipoints.append(NBi)
    Tlpoints.append(NTl)
    Pbpoints.append(NPb)
    Bi209points.append(NBi_209)

    # number of Pb atoms that decay
    Pb_decay = 0
    for i in range(NPb):
        if random.random() < p_Pb:
            Pb_decay += 1
    # add / subtract accordingly
    NPb -= Pb_decay
    NBi_209 += Pb_decay

    # number of Tl atoms that decay
    Tl_decay = 0
    for i in range(NTl):
        if random.random() < p_Tl:
            Tl_decay += 1
    # add / subtract accordingly
    NTl -= Tl_decay
    NPb += Tl_decay

    # number of Bi atoms that decay
    Bi_decay = 0
    for i in range(NBi):
        if random.random() < p_Bi:
            Bi_decay += 1
    # decide which atom the Bi decays into
    Bi_to_Tl = 0
    Bi_to_Pb = 0
    for i in range(Bi_decay):
        if random.random() < 0.9791:
            Bi_to_Pb += 1
        else:
            Bi_to_Tl += 1
    # add / subtract accordingly
    NBi -= Bi_decay
    NTl += Bi_to_Tl
    NPb += Bi_to_Pb

# create plot
plt.plot(tpoints, Bipoints, 'k-', label='Bi 213')
plt.plot(tpoints, Tlpoints, 'r-', label='Tl')
plt.plot(tpoints, Pbpoints, 'g-', label='Pb')
plt.plot(tpoints, Bi209points, 'b-', label='Bi 209')
plt.legend()
plt.title('Radioactive Decay')
plt.xlabel('Time (s)')
plt.ylabel('Number of atoms')
plt.savefig('decaychain.png')
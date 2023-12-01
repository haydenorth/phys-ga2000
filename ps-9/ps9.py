# Graduate Computational Physics
# Problem Set 9
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

# spatial points
L = 1
x = np.linspace(0, L, 1001)
a = L/1000

# time points
h = 10**-4
time = np.arange(0, 0.1, h)

# initial condition: psi(t = 0)
x0 = L/2
sigma = L/10
kappa = 50/L
psi = np.exp(-1* ((x - x0)**2)/(2*sigma**2)) * np.exp(1j * kappa * x)
psi[0] = psi[-1] = 0

# create banded matrix A
N = 1001
a1 = 1 + 1j*(h/(a**2))
a2 = -1j*(h/(2*a**2))
# create diagonals of banded matrix
A_main_diag = np.ones(N, dtype=complex) * a1
# upper diagonal
A_upper = np.ones(N, dtype=complex) * a2
A_upper[0] = 0
# lower diagonal
A_lower = np.ones(N, dtype=complex) * a2
A_lower[-1] = 0
# build matrix
A = np.array([A_upper, A_main_diag, A_lower])

# b1 and 2
b1 = 1 - 1j*(h/(a**2))
b2 = 1j*(h/(2*a**2))

# solve banded matrix problem repeatedly to get psi values
psi_values = []
for t in time:
    # make copies
    psi_values.append(psi.copy())
    psiold = psi.copy()
    
    # calculate vector v
    # psi does not contain the endpoints, so need to append the boundary condition
    psiold = np.concatenate(([0],psi,[0])) 
    v = b1*psiold[1:-1] + b2*(psiold[2:]+psiold[:-2])
    
    # solve banded matrix
    psi = linalg.solve_banded((1,1), A, v)
    psi[0] = psi[-1] = 0

# convert psi values list to array
psi_values = np.array(psi_values, dtype=complex)
# take the real parts
real_parts = np.real(psi_values)

# make a few plots at different times
# make it obvious that the wavefunction moves to the right then spreads out
#for i in np.arange(0, 40):
#    plt.plot(x, real_parts[i], label='Psi(t=' + str(time[i]) + ')')
def plot1(x, real_parts):
    plt.plot(x, real_parts[0], label='Psi(t=' + str(time[0]) + ')')
    plt.plot(x, real_parts[40], label='Psi(t=' + str(time[40]) + ')')
    plt.plot(x, real_parts[110], label='Psi(t=' + str(time[110]) + ')')
    plt.title('Real Wavefunction Psi(t) with hbar/2m = 1')
    plt.xlabel('Position * 1e-8 (m)')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.savefig('ps9_plot1')

def plot2(x, real_parts):
    plt.plot(x, real_parts[999], label='Psi(t=' + str(time[999]) + ')')
    plt.title('Real Wavefunction Psi(t) with hbar/2m = 1')
    plt.xlabel('Position * 1e-8 (m)')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.savefig('ps9_plot4')

plot2(x, real_parts)
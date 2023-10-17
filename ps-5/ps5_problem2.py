# Graduate Computational Physics
# Problem Set 5
# Problem 2
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt

# code to import data from .dat file
# function for checking if read data is a float
def is_float(string):
    """ True if given string is float else False"""
    try:
        return float(string)
    except ValueError:
        return False
# read in data from signal.dat
data = []
with open('signal.dat', 'r') as f:
    d = f.readlines()
    for i in d:
        k = i.rstrip().split('|')
        for i in k:
            if is_float(i):
                data.append(float(i)) 

# clean up data read in from signal.dat
data = np.array(data, dtype='float')
time = data[::2]
signal = data[1::2]
time = time*1e-9

# part a
# plot data
def plot_raw(t, s):
    plt.plot(t, s, '.')
    plt.title('Signal Data')
    plt.xlabel('Time (s)')
    plt.ylabel('Signal')
    plt.savefig('ps5_p2_a')

#plot_raw(time, signal)

# part b
# use SVD to find best third order fit
# create A matrix
A = np.zeros((len(signal), 4))
A[:, 0] = 1.
A[:, 1] = time
A[:, 2] = time**2
A[:, 3] = time**3
# SVD it
(u, w, vt) = np.linalg.svd(A, full_matrices=False)
# invert A
Ainv = vt.transpose().dot(np.diag(1. / w)).dot(u.transpose())
# multiply Ainv*signal
c = Ainv.dot(signal)
# get signal fit
signal_fit = A.dot(c)
'''
plt.plot(time, signal, 'b.', label='Signal')
plt.plot(time, signal_fit, 'r.', label='Cubic Fit')
plt.legend()
plt.title('Signal with Cubic Fit')
plt.xlabel('Time *10^(-9) (s)')
plt.ylabel('Signal')
plt.savefig('ps5_p5_b')
'''

# part c calculate residuals
residuals = signal - signal_fit
'''
plt.plot(time, residuals, '.')
plt.title('Residuals of Cubic Fit')
plt.xlabel('Time *10^(-9) (s)')
plt.ylabel('Residual')
plt.savefig('ps5_p2_c')
'''
# the plot of residuals looks the same as the original plot lol

# part d higher order poly
# create A matrix
A = np.zeros((len(signal), 9))
A[:, 0] = 1.
A[:, 1] = time
A[:, 2] = time**2
A[:, 3] = time**3
A[:, 4] = time**4
A[:, 5] = time**5
A[:, 6] = time**6
A[:, 7] = time**7
A[:, 8] = time**8

# SVD it
(u, w, vt) = np.linalg.svd(A, full_matrices=False)
# invert A
Ainv = vt.transpose().dot(np.diag(1. / w)).dot(u.transpose())
# multiply Ainv*signal
c = Ainv.dot(signal)
# get signal fit
signal_fit = A.dot(c)
# print singular numbers to get condition number 
print(w)
'''
plt.plot(time, signal, 'b.', label='Signal')
plt.plot(time, signal_fit, 'r.', label='Order-8 Fit')
plt.legend()
plt.title('Signal with Polynomial Fit')
plt.xlabel('Time *10^(-9) (s)')
plt.ylabel('Signal')
plt.savefig('ps5_p2_d')
'''
# ok it looks like a polynomial fit isnt going to work :-/
# condition number keeps getting larger as order increases and fit isn't getting better

# part e 
# now try a fourier series type fit

# get linear trend and subtract it off
A = np.zeros((len(time), 2))
A[:, 0] = 1.
A[:, 1] = time
(u, w, vt) = np.linalg.svd(A, full_matrices=False)
ainv = vt.transpose().dot(np.diag(1. / w)).dot(u.transpose())
c = ainv.dot(signal)
ym = A.dot(c) 
# subtract it off
flat_signal = signal-ym
#plt.plot(time, flat_signal, '.')
#plt.show()

# now fit sin/cos to this flat signal
T = 0.5 # half of the time span covered
freq = (2*np.pi) / T
# use 6 integer harmonics of this frequency
A = np.zeros((len(time), 13))
A[:, 0] = 1.
A[:, 1] = np.sin(freq*time)
A[:, 2] = np.cos(freq*time)
A[:, 3] = np.sin(2*freq*time)
A[:, 4] = np.cos(2*freq*time)
A[:, 5] = np.sin(3*freq*time)
A[:, 6] = np.cos(3*freq*time)
A[:, 7] = np.sin(4*freq*time)
A[:, 8] = np.cos(4*freq*time)
A[:, 9] = np.sin(5*freq*time)
A[:, 10] = np.cos(5*freq*time)
A[:, 11] = np.sin(6*freq*time)
A[:, 12] = np.cos(6*freq*time)

(u, w, vt) = np.linalg.svd(A, full_matrices=False)
ainv = vt.transpose().dot(np.diag(1. / w)).dot(u.transpose())
c = ainv.dot(flat_signal)
fourier_fit = A.dot(c)
#print(w)
# plot
'''
plt.plot(time, flat_signal, 'b.', label='Signal')
plt.plot(time, 2*fourier_fit, 'r.', label='Fit')
plt.legend()
plt.title('Signal with Lomb-Scargle Fit')
plt.xlabel('Time *10^(-9) (s)')
plt.ylabel('Signal')
plt.savefig('ps5_p2_e')
'''

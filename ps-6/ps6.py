# Graduate Computational Physics
# Problem Set 6
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy.linalg as linalg
import time

# part a
# read in data from file
hdu_list = fits.open('specgrid.fits')
logwave = hdu_list['LOGWAVE'].data
flux = hdu_list['FLUX'].data

# wavelengths in nm 
wave = (10**logwave) / 10

# function for making the part a plot: plot of the spectra of 4 galaxies
def plot_parta(wave, flux):
    plt.plot(wave, flux[1], label='Galaxy 1')
    plt.plot(wave, flux[2], label='Galaxy 2')
    plt.plot(wave, flux[3], label='Galaxy 3')
    plt.plot(wave, flux[4], label='Galaxy 4')
    plt.legend()
    plt.title('Galaxy Spectra')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (10^-17 erg s^-1 cm^-2 A^-1)')
    plt.savefig('ps6_parta')
#plot_parta(wave, flux)

# part b
# normalize galaxy fluxes
Ngalaxy, Nspectrum = flux.shape
normalized_flux = np.zeros(flux.shape)
normalizations = []

for i in range(Ngalaxy):
    flux_sum = np.sum(flux[i][:])
    norm_coeff = 1/flux_sum
    normalizations.append(norm_coeff)
    normalized_flux[i][:] = norm_coeff*flux[i][:]

# part c
# calculate means and create residuals matrix
i_means = []
residuals = np.zeros(flux.shape)
for i in range(Ngalaxy):
    i_mean = np.average(normalized_flux[i][:])
    residuals[i][:] = normalized_flux[i][:] - i_mean
    i_means.append(i_mean)

# part d
time1 = time.time() # for timing
# calculate covariance matrix C
C = residuals.transpose().dot(residuals)
# calculate eigenvalues and eigenvectors of C
a, v = linalg.eig(C)
time2 = time.time() # timing
print(f'Covariant matrix took {time2-time1} s')
print(a[0:5])

# function to plot eigenvectors
def plot_partd(v):
    plt.plot(v[0], label='Vector 1')
    plt.plot(v[1], label='Vector 2')
    plt.plot(v[2], label='Vector 3')
    plt.plot(v[3], label='Vector 4')
    plt.plot(v[4], label='Vector 5')
    plt.title('Eigenvectors of Covariance Matrix')
    plt.xlabel('Eigenvector index')
    plt.ylabel('Eigenvector value')
    plt.legend()
    plt.savefig('ps6_partd')
#plot_partd(v)

# part e
# find eigenvalues/vectors via svd method
time3 = time.time() # for timing
(u, w, vt) = linalg.svd(residuals)
time4 = time.time() # timing
print(f'SVD took {time4-time3} s')
#w = w**2
# the eigenvalues check out
#print(w[0:5])
# show the vectors are equivalent to what I found before (show first 10)
V = vt.transpose() # eigenvectors are stored in V
for i in range(10):
    plt.plot(v[:,i], V[:,i])
plt.title('Covariance vs SVD Eigenvectors')
plt.xlabel('Covariance Matrix Eigenvectors')
plt.ylabel('SVD Eigenvectors')
plt.savefig('ps6_parte_2.png')

# function to plot the eigenvectors
def plot_parte(v):
    plt.plot(v[0], label='Vector 1')
    plt.plot(v[1], label='Vector 2')
    plt.plot(v[2], label='Vector 3')
    plt.plot(v[3], label='Vector 4')
    plt.plot(v[4], label='Vector 5')
    plt.title('Eigenvectors via SVD')
    plt.xlabel('Eigenvector index')
    plt.ylabel('Eigenvector value')
    plt.legend()
    plt.savefig('ps6_parte')
# take transpose of eigenvector matrix
#plot_parte(V)

# part f, compare the methods by comparing the condition numbers
# C condition number
a = np.sort(np.abs(a))
Ccond = a[-1] / a[0]
print(Ccond)
# R condition number
w = np.sort(w)
Rcond = w[-1] / w[0]
print(Rcond)

# part g
# use eigenvectors to approximate the data: get coefficients of each eigenvector to
# represent the data. Use only the first 5 eigenvectors
eigvec=V[:,:5]
reduced_wavelength_data= np.dot(eigvec.T,residuals.T)
# calculate eigenvector weights
reduced_weights = reduced_wavelength_data.T
# calculate approximation of data
approximate_data = np.dot(eigvec, reduced_wavelength_data).T

# function for plotting data approximized with 5 eigenvectors
def plot_partg(logwave, approximate_data):
    plt.plot(logwave, approximate_data[0,:])
    plt.plot(logwave, approximate_data[1,:])
    plt.plot(logwave, approximate_data[2,:])
    plt.plot(logwave, approximate_data[3,:])
    plt.plot(logwave, approximate_data[4,:])
    plt.title('Approximate spectra with 5 eigenvectors')
    plt.xlabel('Log(wavelength) (A)')
    plt.ylabel('Normalized Flux (10^-17 erg s^-1 cm^-2 A^-1)')
    plt.savefig('ps6_partg_2')
#plot_partg(logwave, approximate_data)

# part h
# plot c1 and c2 vs c0
def plot_parth(reduced_weights):
    plt.plot(reduced_weights[:,0], reduced_weights[:,1], 'o', label='c1')
    plt.plot(reduced_weights[:,0], reduced_weights[:,2], 'o', label='c2')
    plt.title('Eigenvector weights')
    plt.xlabel('c0')
    plt.ylabel('c1 and c2')
    plt.legend()
    plt.savefig('ps6_parth.png')
#plot_parth(reduced_weights)

# part i
# calculate squared fractional residual for varying number of eigenvectors N
total_sq_resids = []
for i in range(20):
    eigvec=V[:,:i]
    reduced_wavelength_data = np.dot(eigvec.T,residuals.T)
    approximate_data = np.dot(eigvec, reduced_wavelength_data).T
    # above is same as before
    # no calculate the residual squared
    squared_residual = (approximate_data - normalized_flux)**2
    total_sq_resids.append(np.sum(squared_residual))
# plot total residual as a function of N
def plot_parti(total_sq_resids):
    plt.plot(range(20), total_sq_resids)
    plt.title('Total Squared Residuals')
    plt.xlabel('Number of eigenvectors, N')
    plt.ylabel('Total Squared Residual')
    plt.savefig('ps6_parti')
#plot_parti(total_sq_resids)


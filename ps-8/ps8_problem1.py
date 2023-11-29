# Graduate Computational Physics
# Problem Set 8
# Problem 1
# By Hayden Orth

import numpy as np
import matplotlib.pyplot as plt
from scipy import fft

# read data into arrays from txt files
trumpet = np.loadtxt('trumpet.txt')
piano = np.loadtxt('piano.txt')

# function for plotting
def plot_waveform(waveform):
    # plot waveform
    plt.plot(waveform)
    plt.title('Piano Waveform')
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.savefig('ps8_p1_a_piano.png')

# function for plotting fft
def plot_fft(waveform):
    # take fourier transform
    c = fft.fft(waveform)
    # take magnitude of coefficients
    mag = abs(c)**2
    # take only first 10000 coefficients
    #mag = mag[:10000]

    # calculate frequency resolution
    # number of samples N
    N = len(waveform)
    # Nyquist frequency = 1/2*(spacing)
    fNy = 1/(2*(1/44100))
    # frequency resolution = 1/N*(spacing)
    freq_res = 1/(N*(1/44100))

    # calculate frequency values
    freq = np.arange(0, 2*fNy, freq_res)

    # plot spectrum (first 10000 coefficients)
    plt.plot(freq[:10000], mag[:10000])
    plt.title('Piano Fourier Spectrum')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')
    #plt.vlines(523, 0, 10e15, colors='k')
    # the main peak is at 523 Hz which is C5
    plt.savefig('piano_fourier_spectrum')

#plot_waveform(piano)
plot_fft(piano)
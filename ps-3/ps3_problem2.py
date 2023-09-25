# Graduate Computational Physics
# Problem Set 3
# Problem 2
# By Hayden Orth

import numpy as np
import time
import matplotlib.pyplot as plt

# function for multiplying two NxN matrices, returns resultant matrix
def mult_matrix(A, B, N):
    C = np.zeros([N,N], float)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i,j] += A[i,k]*B[k,j]
    return C

# gather data: computation time as a function of matrix size N
Ns = np.arange(50, 550, 50)
times1 = []
times2 = []
# test using mult_matrix function and dot() method
for N in Ns:
    A = np.ones([N,N])*24
    B = np.ones([N,N])*420
    ti = time.time()
    mult_matrix(A, B, N)
    tf = time.time()
    times1.append(tf-ti)
    ti = time.time()
    np.dot(A, B)
    tf = time.time()
    times2.append(tf-ti)
    # to keep track of where we're at:
    print(N)

# create plot
times1 = np.array(times1)
times2 = np.array(times2)
plt.plot(Ns, times1, 'r-', label='For loops')
plt.plot(Ns, times2, 'b-', label='NumPy dot() method')
plt.title('Matrix Multiplication')
plt.xlabel('Matrix size NxN')
plt.ylabel('Computation time (s)')
plt.legend()
plt.savefig('matrixmult.png')

    


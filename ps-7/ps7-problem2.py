# Graduate Computational Physics
# Problem Set 7
# Problem 2
# By Hayden Orth

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import optimize

# read in data with pandas
data = pd.read_csv('survey.csv')
# convert to numpy arrays
xs = data['age'].to_numpy()
ys = data['recognized_it'].to_numpy()
x_sort = np.argsort(xs)
age = xs[x_sort]
recognize = ys[x_sort]

# function for probability
def p(x, beta0, beta1): return 1/(1+np.exp(-(beta0+beta1*x)))

# function for log likelihood
def log_likelihood(beta, xs, ys):
    beta_0 = beta[0]
    beta_1 = beta[1]
    epsilon = 1e-16
    l_list = [ys[i]*np.log(p(xs[i], beta_0, beta_1)/(1-p(xs[i], beta_0, beta_1)+epsilon)) 
              + np.log(1-p(xs[i], beta_0, beta_1)+epsilon) for i in range(len(xs))]
    ll = np.sum(np.array(l_list), axis = -1)
    return -ll # return log likelihood

# Now maximize log likelihood
pstart = [1,42]

# Covariance matrix of parameters
def Covariance(hess_inv, resVariance):
    return hess_inv * resVariance

# Error of parameters
def error(hess_inv, resVariance):
    covariance = Covariance(hess_inv, resVariance)
    return np.sqrt( np.diag( covariance ))

# maximize log likelihood
b = np.linspace(-4, 4, 100)
result = optimize.minimize(lambda b, age, recognize: log_likelihood(b, age, recognize), [0,0], args=(age, recognize))
hess_inv = result.hess_inv # inverse of hessian matrix
var = result.fun/(len(recognize)-len(pstart)) 
dFit = error( hess_inv,  var)
print('Optimal parameters and error:\n\tp: ' , result.x, '\n\tdp: ', dFit)
print('Covariance matrix of optimal parameters:\n\tC: ' , Covariance( hess_inv,  var))

# plot data and optimized logistic function
plt.plot(age, p(age, result.x[0], result.x[1]), label='Logistic Function')
plt.plot(age, recognize, 'o', label='Data')
plt.title('Probability of hearing the phrase')
plt.xlabel('Age (years)')
plt.ylabel('Probabilty')
plt.legend()
plt.savefig('ps7.png')

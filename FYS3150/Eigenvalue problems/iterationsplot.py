import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit

data = np.loadtxt('iterations.txt')
N = data[:, 0]
ite = data[:, 1]

ns = np.linspace(np.min(N), np.max(N), 100)
y = ns**2.03

plt.plot(N, ite, 'o', label='Data points')
plt.plot(ns, y, label='$f(x) = N^{2.03}$')
plt.xlabel('Size of matrix N')
plt.ylabel('Iterations')
plt.legend()
plt.savefig('Iterationsplot.pdf')

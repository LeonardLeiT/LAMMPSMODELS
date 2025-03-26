#coding = utf - 8


import numpy as np 
import matplotlib.pyplot as plt 

def potential(epsilon, sigma, r, alpha = 2):
    return np.where(r <= sigma, epsilon / alpha * (1 - r / sigma)**alpha, 0)

def force(epsilon, sigma, r, alpha = 2):
    return np.where(r <= sigma, epsilon / sigma * (1 - r / sigma) ** (alpha - 1), 0)

r = np.linspace(0.5, 1.5, 100)

plt.figure(figsize = (10, 5))
plt.subplot(121)
data = np.loadtxt('table11.dat', skiprows = 5)
plt.plot(data[:, 1], data[:, 2], 'o', label = 'Pair_Write')
plt.plot(r, potential(1.0, 1.0, r), lw = 2, label = 'function')
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$V(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
plt.legend(fontsize = 15)

plt.subplot(122)
plt.plot(data[:, 1], data[:, 3], 'o')
plt.plot(r, force(1.0, 1.0, r), lw = 2)
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$F(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)

plt.tight_layout()
plt.savefig('compare.png')
plt.show()
plt.close()
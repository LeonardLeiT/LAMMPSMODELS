#coding = utf - 8


import numpy as np 
import matplotlib.pyplot as plt 

def potential(epsilon, sigma, r, alpha = 2):
    return np.where(r <= sigma, epsilon / alpha * (1 - r / sigma)**alpha, 0)

def force(epsilon, sigma, r, alpha = 2):
    return np.where(r <= sigma, epsilon / sigma * (1 - r / sigma) ** (alpha - 1), 0)

r = np.linspace(0.2, 1.5, 100)
plt.plot(r, force(1.0, 1.0, r), label = 'function', zorder = 10, lw = 3)

distance = []
force = []
for n in range(50):
    n += 1
    data = np.loadtxt('dump.n%d.atom'%n, skiprows = 10)[2:]
    distance.append(np.linalg.norm(data[:3]))
    force.append(np.linalg.norm(data[3:]))

plt.plot(distance, force, 'o', label = 'lammps-two atoms')


plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$F(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
plt.legend(fontsize = 16)

plt.tight_layout()
plt.savefig('compare.png')
plt.show()
plt.close()

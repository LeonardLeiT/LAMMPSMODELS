#coding = utf - 8


import numpy as np 
import matplotlib.pyplot as plt 

def funct_form(A, epsilon, sigma, r, n):
    E = A*epsilon*np.power(sigma/r, n)
    F = A*epsilon*n/r*np.power(sigma/r, n)
    return E, F

def energy_force(A, epsilon, sigma, r, n, r_cut=3.0):
    r_cut *= sigma
    E, F = funct_form(A, epsilon, sigma, r, n)
    Ecut, Fcut = funct_form(A, epsilon, sigma, r_cut, n)
    energy = E - Ecut + (r - r_cut)*Fcut
    force  = F - Fcut
    condition = (r > r_cut)
    energy[condition] = 0
    force[condition] = 0
    return energy, force

plt.figure(figsize=(6, 6))
plt.subplot(211)
r = np.linspace(0.2, 2.0, 180)
epsilon = 1.0
sigma = 1.4
A = 1.0
powern = 36
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut=3.0)
plt.plot(r, force, label = 'function', zorder = 10, lw = 3)

distance = []
forceMD = []
energyMD = []
for n in range(200, 600):
    n += 1
    data = np.loadtxt('dump.n%d.atom'%n, skiprows = 21)[2:]
    distance.append(np.linalg.norm(data[:3]))
    forceMD.append(np.linalg.norm(data[3:-1]))
    energyMD.append(data[-1]*2)


plt.plot(distance, forceMD, 'o', label = 'lammps-two atoms')

plt.plot([0.5, 4], [0, 0], '--')

plt.ylim(-1.2, 3.0)
plt.xlim(0.8, 2)
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$F(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
plt.title(r'$\epsilon=%s, \sigma=%s, A=%s, powern=%s$'%(epsilon, sigma, A, powern), size=20)
plt.legend(fontsize = 13, loc='upper left')

plt.subplot(212)
plt.plot(r, energy, label = 'function', zorder = 10, lw = 3)
plt.plot(distance, energyMD, 'o', label = 'lammps-two atoms')
plt.plot([0.5, 4], [0, 0], '--')
plt.ylim(-1.2, 3.0)
plt.xlim(0.8, 2)
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$E(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
#plt.legend(fontsize = 16)

plt.tight_layout()
plt.savefig('compare_force_22.png')
plt.show()
plt.close()

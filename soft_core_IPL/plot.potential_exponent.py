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
r = np.linspace(0.5, 1.5, 200)
epsilon = 1.0
sigma = 1.0
A = 1.0
r_cut = 1.3 * sigma

powern = 36
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, force, label = r'$n=%d$'%powern, lw = 2)

powern = 27
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, force, label = r'$n=%d$'%powern, lw = 2)

powern = 21
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, force, label = r'$n=%d$'%powern, lw = 2)

powern = 12
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, force, label = r'$n=%d$'%powern, lw = 2)

plt.plot([0.5, 4], [0, 0], '--')

plt.ylim(-1.2, 3.0)
plt.xlim(0.8, 1.5)
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$F(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
plt.title(r'$\epsilon=%s, \sigma=%s, A=%s, r_{cut}=%s$'%(epsilon, sigma, A, r_cut), size=20)
plt.legend(fontsize = 13, loc='upper left')

plt.subplot(212)
powern = 36
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, energy, label = r'$n=%d$'%powern, lw = 2)

powern = 27
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, energy, label = r'$n=%d$'%powern, lw = 2)

powern = 21
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, energy, label = r'$n=%d$'%powern, lw = 2)

powern = 12
energy, force = energy_force(A, epsilon, sigma, r, powern, r_cut)
plt.plot(r, energy, label = r'$n=%d$'%powern, lw = 2)


plt.plot([0.5, 4], [0, 0], '--')
plt.ylim(-1.2, 3.0)
plt.xlim(0.8, 1.5)
plt.xlabel(r'$r$', size = 20)
plt.ylabel(r'$E(r)$', size = 20)
plt.xticks(size = 18)
plt.yticks(size = 18)
#plt.legend(fontsize = 16)

plt.tight_layout()
plt.savefig('potential_exponent.png')
plt.show()
plt.close()

import numpy as np 
import matplotlib.pyplot as plt 

def potential_force(r):
    m = 16
    A = 5.82
    c = 1.1
    a = 1.87
    B = 1.28
    d = 0.27
    b = 1.94

    prefactor = A * np.exp(c/(r-a))
    V1 = prefactor * (r**(-m) - B)
    dV1 = prefactor * (-m * r**(-m-1)) + V1*(-c/(r-a)**2)
    V1 = np.where(r<a, V1, 0)
    dV1 = np.where(r<a, dV1, 0)

    V2 = B * np.exp(d/(r-b))
    dV2 = V2 * (-d/(r-b)**2)
    V2 = np.where(r<b, V2, 0)
    dV2 = np.where(r<b, dV2, 0)

    potential = V1 + V2
    force = dV1 + dV2
    return potential, -force

r = np.linspace(0.5, 3, 1000)
potential, force = potential_force(r)
plt.figure(figsize=(12, 5))

plt.subplot(121)
plt.plot(r, potential, "-o")
plt.ylim(-1, 3)
plt.xlim(0.6, 2.2)
plt.xlabel("r")
plt.ylabel(r"V(r)")

plt.subplot(122)
plt.plot(r, force, "-o")
plt.ylim(-5, 5)
plt.xlim(0.6, 2.2)
plt.xlabel("r")
plt.ylabel(r"F(r)")

plt.show()
plt.close()
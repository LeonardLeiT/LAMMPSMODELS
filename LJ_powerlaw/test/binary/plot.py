#coding = utf-8

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

filename = 'dynamics.linear.T0.002.n1.dat'
data = pd.read_csv(filename, sep = '\s+')
plt.plot(data['t'], data['msd'], '-o', label = 'lammps')

filename = '../tonghua/MSD_Binary_4096_0.91_2.0d-3_2D_T0.002.dat'
data = np.loadtxt(filename)
plt.plot(data[:, 0], data[:, 1], label = 'tonghua')


plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$t$', size = 20)
plt.ylabel(r'$MSD$', size = 20)
plt.legend(fontsize = 18)
plt.tight_layout()
plt.savefig('msd_compare.png')
plt.show()
plt.close()
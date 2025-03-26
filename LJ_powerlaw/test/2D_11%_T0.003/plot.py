#coding = utf-8

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

filename = 'dynamics.log.T0.0030.n1.dat'
data = pd.read_csv(filename, sep = '\s+')
plt.plot(data['t'], data['ISF'], '-o', label = 'direct')

filename = 'dynamics.cage.log.T0.0030.n1.dat'
data = pd.read_csv(filename, sep = '\s+')
plt.plot(data['t'], data['ISF'], '-s', label = 'cage relative')

plt.title('2D_11%_T0.003', size = 20)

plt.xscale('log')
plt.xlabel(r'$t$', size = 20)
plt.ylabel(r'$F_s(q,t)$', size = 20)
plt.legend(fontsize = 18)
plt.tight_layout()
plt.savefig('ISF.png')
plt.show()
plt.close()
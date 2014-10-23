#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

# Run the executable integrate_bernus.out
#os.system("./integrate_bernus.out")

# Read plain text data
data = np.loadtxt("out.txt")
t    = data[:,0]
i_na = data[:,1]
i_ca = data[:,2]
i_to = data[:,3]
i_k  = data[:,4]
i_k1 = data[:,5]
i_b_ca = data[:,6]
i_b_na = data[:,7]
i_na_k = data[:,8]
i_na_ca = data[:,9]
v0 = data[:,10]
Iion = data[:,11]

# And plot it.
plt.figure(figsize=(8,8))
plt.plot(t, i_na, label='i_na', linewidth=2.0)
plt.plot(t, i_ca, label='i_ca', linewidth=2.0)
plt.plot(t, i_to, label='i_to', linewidth=2.0)
plt.plot(t, i_k,  label='i_k', linewidth=2.0)
plt.plot(t, i_k1, label='i_k1', linewidth=2.0)
plt.plot(t, i_b_ca, label='i_b_ca', linewidth=2.0)
plt.plot(t, i_b_na, '--', label='i_b_na', linewidth=2.0)
plt.plot(t, i_na_k, '--', label='i_na_k', linewidth=2.0)
plt.plot(t, i_na_ca, '--', label='i_na_ca', linewidth=2.0)
plt.plot(t, v0, '--', label='V', linewidth=2.0)
plt.plot(t, Iion, '--', label='Iion', linewidth=2.0)

plt.legend()
plt.xlim(t[0], t[np.size(t)-1])

plt.xlabel('Time [ms]', fontsize=14)
#plt.ylabel('Potential [mV]', fontsize=14)
plt.show()

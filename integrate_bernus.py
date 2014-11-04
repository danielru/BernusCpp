#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

# Run the executable integrate_bernus.out
os.system("./integrate_bernus.out")

# Read plain text data
data = np.loadtxt("bernus.txt")
t    = data[:,0]
V0   = data[:,1]
m    = data[:,2]
v    = data[:,3]
f    = data[:,4]
to   = data[:,5]
x    = data[:,6]
Iion = data[:,7]

# And plot it.
plt.figure(figsize=(8,8))
plt.plot(t, V0, label='V', linewidth=2.0)
plt.xlim(t[0], t[np.size(t)-1])
plt.ylim((-100, 60))
plt.xlabel('Time [ms]', fontsize=14)
plt.ylabel('Potential [mV]', fontsize=14)
plt.show()
#plt.savefig('potential.eps')

# And plot it.
plt.figure(figsize=(8,8))
plt.plot(t, Iion, label='Iion', linewidth=2.0)
plt.xlim(t[0], t[np.size(t)-1])
plt.xlabel('Time [ms]', fontsize=14)
plt.ylabel('Current [XX]', fontsize=14)
plt.show()
#plt.savefig('potential.eps')

plt.figure(figsize=(8,8))
plt.plot(t, m, label='m', linewidth=2.0)
plt.plot(t, v, label='v', linewidth=2.0)
plt.plot(t, f, label='f', linewidth=2.0)
plt.plot(t, to, label='to', linewidth=2.0)
plt.plot(t, x, label='x', linewidth=2.0)
plt.legend()
plt.xlim(t[0], t[np.size(t)-1])
plt.xlabel('Time [ms]', fontsize=14)
plt.ylabel('Gate variable', fontsize=14)
plt.show()
#plt.savefig('gates.eps')
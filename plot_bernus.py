#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

os.system("./probe_bernus.out")

#np.loadtxt("vanderpol.txt")
data = np.loadtxt("bernus.txt")
t  = data[:,0]
V0 = data[:,1]
m  = data[:,2]
v  = data[:,3]
f  = data[:,4]
to = data[:,5]
x  = data[:,6]
Iion = data[:,7]

plt.figure(figsize=(8,8))
plt.plot(t, V0, label='V')
plt.plot(t, Iion, label='Ion')
plt.legend()
plt.xlim(t[0], t[np.size(t)-1])
plt.ylim((-100, 60))
plt.show()

plt.figure(figsize=(8,8))
plt.plot(t, m, label='m')
plt.plot(t, v, label='v')
plt.plot(t, f, label='f')
plt.plot(t, to, label='to')
plt.plot(t, x, label='x')
plt.legend()
plt.xlim(t[0], t[np.size(t)-1])
plt.xlabel("Time (ms)", fontsize=14)
plt.ylabel('Gate variable', fontsize=14)
plt.show()

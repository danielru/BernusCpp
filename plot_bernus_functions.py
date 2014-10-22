#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

os.system("./probe_bernus_functions.out")

#np.loadtxt("vanderpol.txt")
data = np.loadtxt("bernus_functions.txt")
V     = data[:,0]
v_inf = data[:,1]
tau_v = data[:,2]
x_inf = data[:,3]
tau_x = data[:,4]

plt.figure(figsize=(8,8))
plt.plot(V, v_inf)
plt.xlim(V[0], V[np.size(V)-1])
plt.xlabel("membrane potential (mV)", fontsize=14)
plt.ylabel("$v_{\infty}$", fontsize=14)
plt.show()

plt.figure(figsize=(8,8))
plt.plot(V, tau_v)
plt.xlim(V[0], V[np.size(V)-1])
plt.xlabel("membrane potential (mV)", fontsize=14)
plt.ylabel("$\tau_{v}$", fontsize=14)
plt.show()

plt.figure(figsize=(8,8))
plt.plot(V, x_inf)
plt.xlim(V[0], V[np.size(V)-1])
plt.xlabel("membrane potential (mV)", fontsize=14)
plt.ylabel("$x_{\infty}$", fontsize=14)
plt.show()

plt.figure(figsize=(8,8))
plt.plot(V, tau_x)
plt.xlim(V[0], V[np.size(V)-1])
plt.xlabel("membrane potential (mV)", fontsize=14)
plt.ylabel("$\tau_{x}$", fontsize=14)
plt.show()
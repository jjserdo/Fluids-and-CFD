# Imported from CFD Project
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import norm 
##################
start_time = time.time()
s = 5
g = 500
num = 20

X = np.linspace(-s,s,g)
Y = np.linspace(-s,s,g)
x,y = np.meshgrid(X,Y)
r = np.arctan2(y,x)
t = np.sqrt(x**2+y**2)
z = x + 1j*y
''' Input Function Here '''
U = 1
a = np.pi/3
c = np.exp(-1j*a)
Q = 1
G = 1
M = 1
# Success
# z, z**2, c*z, Q/(2*np.pi)*np.log(z), -1j*G/(2*np.pi)*np.log(z)
# U*z**(1.5)

# Weird
# np.exp(z), U*z**(0.25), U*z**(0.75)
Fz = U*c
phi = Fz.real
psi = Fz.imag

fig, ax = plt.subplots(figsize=(5,5))
ax.grid()
ax.set_title('Complex Plots')
ax.set_xlabel('Real')
ax.set_ylabel('Imaginary')
ax.set_xlim([-s,s])
ax.set_ylim([-s,s])
mean = np.mean(phi)
std = np.std(phi)
levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
ax.contour(X, Y, phi, levels=levels, linestyles='dashed',alpha=0.5)
mean = np.mean(phi)
std = np.std(phi)
levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
ax.contour(X, Y, psi, levels=levels)

print("--- %10s seconds ---" % np.round((time.time() - start_time),4))

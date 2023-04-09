# Imported from CFD Project
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import norm 
##################
start_time = time.time()
s = 50
g = 500
num = 20

X = np.linspace(-s,s,g)
Y = np.linspace(-s,s,g)
x,y = np.meshgrid(X,Y)
r = np.arctan2(y,x)
t = np.sqrt(x**2+y**2)
z = x + 1j*y
''' Input Function Here '''
func = 'lmao'
alf = np.radians(90)
a = 10
l = z/2 + np.sign(z.real) * np.sqrt((z/2)**2-a**2) 
Fz = np.exp(-1j*alf)*l + a**2/(np.exp(-1j*alf)*l) 
phi = Fz.real
psi = Fz.imag

fig, ax = plt.subplots(figsize=(5,5))
ax.grid()
ax.set_title('Complex Plots')
ax.set_xlabel('Real')
ax.set_ylabel('Imaginary')
ax.set_xlim([-s,s])
ax.set_ylim([-s,s])
if func == 'lmao':
    mean = np.mean(phi)
    std = np.std(phi)
    levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
    ax.contour(x, y, phi, levels=levels, linestyles='dashed',alpha=0.5)
    mean = np.mean(psi)
    std = np.std(psi)
    levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
    ax.contour(x, y, psi, levels=levels)
else:
    ax.contour(x, y, phi, num, linestyles='dashed',alpha=0.5)
    ax.contour(x, y, psi, num)

print("--- %10s seconds ---" % np.round((time.time() - start_time),4))

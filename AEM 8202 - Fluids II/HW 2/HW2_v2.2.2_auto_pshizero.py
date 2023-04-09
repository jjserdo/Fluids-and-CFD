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
U = 2/3
a = np.pi/3
c = U*np.exp(-1j*a)
Q = 5
G = 1
M = 1
A = 1
sep = 2
# Success
# z, z**2, c*z, Q/(2*np.pi)*np.log(z), -1j*G/(2*np.pi)*np.log(z)
# U*z**(1.5), M/z, U*z + M/z, U*(z+A**2/z)

# Weird
# np.exp(z), U*z**(0.25), U*z**(0.75)
func = 'lmao'
plotzero = 0
plotstag = True
# Fz = -1j*G/(2*np.pi)*np.log(z**3-1) + U*z**(1.5) # Number 1
# Fz = z + Q/(2*np.pi)*np.log((z+sep)/(z-sep))     # Number 3
# Fz = 
#Fz =  1j*M/(c*(z-1j)) + z**2                     # Number 2
#Fz = U*z + (Q*ep/np.pi) / z
Fz = U*z + Q/(2*np.pi)*np.log((z+sep)/(z-sep))
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
if plotzero == True:
    zero = np.argwhere(np.around(psi,1) == 0)
    zeros = np.zeros((2,zero.shape[0]))
    for i in range(zero.shape[0]):
        ax.plot(x[zero[i,0],zero[i,1]],y[zero[i,0],zero[i,1]],'b*')

print("--- %10s seconds ---" % np.round((time.time() - start_time),4))

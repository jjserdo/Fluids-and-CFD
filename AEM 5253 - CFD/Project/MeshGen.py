# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 12:34:31 2022

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt

#### #### Initializing Values of the Structured Grid #### ####

# Parameter Values of Zeta
    # angle
Z0 = 0
Z1 = 2*np.pi
Znum = 52
# Znum = 50

# Paramater Values of Eta
    # axes diameter
N0 = 0.255 # {set on Nov 30} 
N1 = 3.72  # {set on Dec 1 } 
Nnum = 52

# Create Zeta and Eta Values  then Create Grid
''' #Do not Edit '''
Z = np.linspace(Z0, Z1, Znum+1)
N = np.linspace(N0, N1, Nnum+1)
Zgrid, Ngrid = np.meshgrid(Z, N)

# Draws Structured Grid
fig, ax = plt.subplots()
# Draws vertical lines
for i in np.arange(len(Z)):
    ax.plot(Zgrid.T[i], Ngrid.T[i], c='blue')
# Draws horizontal lines
for j in np.arange(len(N)):
    ax.plot(Zgrid[j], Ngrid[j], c='blue') 
ax.set_title('Unstructured Grid')
ax.set_xlim([Z0,Z1])
ax.set_ylim([N0,N1])
plt.show()

#### #### Mapping Structured Grid to Unstructured Grid #### ####
xgrid = np.zeros((len(Z), len(N)))
ygrid = np.zeros((len(Z), len(N)))

# Mapping
for i in np.arange(len(Z)):
    for j in np.arange(len(N)):
        xgrid[i][j] = np.cosh(N[j]) * np.cos(Z[i])
        ygrid[i][j] = np.sinh(N[j]) * np.sin(Z[i])
        
#### Draws Unstructured Grid
fig, ax = plt.subplots()
# Draws the Ellipse itself
for i in np.arange(len(N)):
    ax.plot(xgrid.T[i], ygrid.T[i], c='green')
for j in np.arange(len(Z)):
    ax.plot(xgrid[j], ygrid[j], c='green') 
# Radial
ax.plot(xgrid[int(np.floor(len(Z)/4))][-2], ygrid[int(np.floor(len(Z)/4))][-2], c='blue', marker=".", markersize=10) 
ax.plot(xgrid[int(np.floor(len(Z)/4))][-1], ygrid[int(np.floor(len(Z)/4))][-1], c='blue', marker=".", markersize=10) 
# Angular
ax.plot(xgrid[0][-1], ygrid[0][-1], c='red', marker=".", markersize=10) 
ax.plot(xgrid[1][-1], ygrid[1][-1], c='red', marker=".", markersize=10) 

""" # Dec 1 Addition for Zooming in """
ax.set_title('Unstructured Grid')
#ax.set_xlim([-2,2])
#ax.set_ylim([-2,2])

plt.show()

""" #November 30 """
#### Compute Jacobian ####
J = np.zeros((len(Z), len(N)))

for i in np.arange(len(Z)):
    for j in np.arange(len(N)):
        J[i][j] = np.cosh(N[j])**2 * np.sin(Z[i])**2 + \
                  np.sinh(N[j])**2 * np.cos(Z[i])**2
                  
Zx = np.zeros((len(Z), len(N)))
Nx = np.zeros((len(Z), len(N)))
Zy = np.zeros((len(Z), len(N)))
Ny = np.zeros((len(Z), len(N)))

for i in np.arange(len(Z)):
    for j in np.arange(len(N)):
        Zx[i][j] =  np.cosh(N[j]) * np.sin(Z[i]) / J[i][j]
        Nx[i][j] = -np.sinh(N[j]) * np.cos(Z[i]) / J[i][j]
        Zy[i][j] = -np.sinh(N[j]) * np.cos(Z[i]) / J[i][j]
        Ny[i][j] = -np.cosh(N[j]) * np.sin(Z[i]) / J[i][j]

from matplotlib import cm
'''
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(xgrid,ygrid,J,cmap=cm.coolwarm,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
'''
# December 1 Modification
fig = plt.figure()
ax = plt.axes(projection='3d')
surf = ax.plot_surface(xgrid,ygrid,J,cmap=cm.coolwarm,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.25, aspect=10, fraction=0.1, orientation='horizontal')
ax.view_init(90, 90)
#ax.set_xlim([-2,2])
#ax.set_ylim([-2,2])
plt.show()

# December 3 Addition
plt.contourf(xgrid, ygrid, J)
plt.xlim([-1.25,1.25])
plt.ylim([-1.25,1.25])
plt.show()

""" # December 1 """
# distance is increased to 20 times: N1 = 3.72
# added axes limits

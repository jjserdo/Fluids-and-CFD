# -*- coding: utf-8 -*-
"""
Created on Mon May  1 22:54:20 2023

@author: jjser
Justine John A. Serdoncillo
Spring 2023
AEM 8202 Fluids II
Final Exam
"""

from potentialFlow_airfoil import airfoil
import numpy as np
import matplotlib.pyplot as plt

# Default Values
c = 1 # chord length 
rho = 1 # density
nu = 0.001 # kinematic viscosity

# %% Problem 1a)

# Initialize airfoils
eight = airfoil(c, 8)
four = airfoil(c, 4)
two = airfoil(c, 2)

# Display m values
print(f"m for 8% thickness airfoil is {eight.m}")
print(f"m for 4% thickness airfoil is {four.m}")
print(f"m for 2% thickness airfoil is {two.m}")

# %% Problem 1b)

# Plotting streamlines
eight.plotStreamlines()

# %% Problem 1c)

# Plotting tangential velocities
eight.plotTangentialVelocities()
four.plotTangentialVelocities()
two.plotTangentialVelocities()

# %% Problem 1d)

eight.rotateFlow(4)
eight.plotStreamlines()

# %% Problem 1e)

lift = eight.getLift_Blasius(rho)
print(f"Lift using Blasius is {lift}")
lift = eight.getLift_Integrate(rho)
print(f"Lift using Numerical Integration is {lift}")

# %% Problem 1f)

# Return back to part c conditions
eight.rotateFlow(0)

eight.doThwaites(nu)
four.doThwaites(nu)
two.doThwaites(nu)

''' ... it would separate/not separate '''
print("For 8% thickness airfoil, it seems like it would  ")
print("For 4% thickness airfoil, it seems like it would  ")
print("For 2% thickness airfoil, it seems like it would  ")

# %% Problem 1g)
theta = np.zeros(2,4)
delta = np.zeros(2,4)

for ii in range(4): # 1, 2, 3, 4
    four.rotateFlow(ii+1)
    t, d = four.doThwaites()
    theta[0, ii] = t
    delta[0, ii] = d
        
for ii in range(4): # 1, 2, 3, 4
    eight.rotateFlow(ii+1)
    t, d = eight.doThwaites()
    theta[1, ii] = t
    delta[1, ii] = d
        
fig, ax = plt.subplots()
ax.plot(range(4), theta[0], '-o', label='4')
ax.plot(range(4), theta[1], '-*', label='8')
ax.legend()

print("Thickness is maybe not too good tbh")






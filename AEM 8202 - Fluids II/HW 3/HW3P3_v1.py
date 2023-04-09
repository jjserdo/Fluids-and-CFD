# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 21:44:28 2023

@author: jjser
Justine John A. Serdoncillo 
AEM 8202 Fluids II
Homework 3 Problem 3 v1
    - Solve the Non-Newtonian Similarity ODE
"""
import numpy as np
import matplotlib.pyplot as plt

# %%
"""
    Personal Created Functions
    - g1, g2 and g3 represents 3 linear ODEs
    - rk4 is the Runge-Kutta 4th Order Method
    - shoot is the shooting method needed
"""
def g1(f, df, ddf, n):
    """Returns the value of the first ODE"""
    return df

def g2(f, df, ddf, n):
    """Returns the value of the second ODE"""
    return ddf

def g3(f, df, ddf, n):
    """Returns the value of the third ODE"""
    return - f * ddf / (ddf**(1/2))

def rk4(ni, nf, f0, df0, ddf0, dn):
    N = int((nf - ni) / dn)
    n = np.linspace(ni, nf, N+1)
    F = np.zeros(len(n))
    dF = np.zeros(len(n))
    ddF = np.zeros(len(n))
    F[0] = f0
    dF[0] = df0
    ddF[0] = ddf0
    for i in range(len(n)-1):
        k01 = g1(F[i], dF[i], ddF[i], n[i])
        k02 = g2(F[i], dF[i], ddF[i], n[i])
        k03 = g3(F[i], dF[i], ddF[i], n[i])
        
        k11 = g1(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, n[i]+dn/2)
        k12 = g2(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, n[i]+dn/2)
        k13 = g3(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, n[i]+dn/2)
        
        k21 = g1(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, n[i]+dn/2)
        k22 = g2(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, n[i]+dn/2)
        k23 = g3(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, n[i]+dn/2)
        
        k31 = g1(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, n[i]+dn)
        k32 = g2(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, n[i]+dn)
        k33 = g3(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, n[i]+dn)
        
        F[i+1] = F[i] + dn/6 * (k01 + 2*k11 + 2*k21 + k31)
        dF[i+1] = dF[i] + dn/6 * (k02 + 2*k12 + 2*k22 + k32)
        ddF[i+1] = ddF[i] + dn/6 * (k03 + 2*k13 + 2*k23 + k33)
    return n, F, dF, ddF


# %% 
"""
   Non-Newtonian Boundary Layer Profile
   - dddf * ddf^(1/2) + f * ddf = 0
"""
# Define problem-specific variables
f0 = 0        # initial value of the dependent variable f
df0 = 0       # initial value of the first derivative of f
ddf0 = 0.617  # initial value of the second derivative of f
ni = 0.0      # initial value of the independent variable n
nf = 5.0      # final value of the independent variable n
dn = 0.01     # step size

# Solve the Falkner-Skan equations using the Runge-Kutta method
n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn)

# Create a plot to display the solution
fig, ax = plt.subplots(figsize=(7,5))

# Set the title and axis labels for the plot
ax.set_title('Comparison betweeen Newtonian and Non-Newtonian \n Boundary Layers')
ax.set_ylabel('$ \\eta $')
ax.set_xlabel('$ f\'(\\eta) $')

# Set the x and y limits for the plot
ax.set_ylim([0,nf])
#ax.set_ylim([0,2.0])

# Plot a horizontal dashed line at y=1 to show the target value
ax.plot(np.ones(100), np.linspace(0,100,100,True), '--')

# Plot the solution for f'(eta)
ax.plot(dF, n, label='Non-Newtonian')

# Add gridlines to the plot
ax.grid()

# %% 
"""
    Blasius Newtonian Solution Taken from Problem 2
"""
# Define personal created functions
def g1(f, df, ddf, a, n):
    """Returns the value of the first ODE"""
    return df

def g2(f, df, ddf, a, n):
    """Returns the value of the second ODE"""
    return ddf

def g3(f, df, ddf, a, n):
    """Returns the value of the third ODE"""
    return -0.5 * (a+1) * f * ddf - a * (1 - df**2)

def rk4(ni, nf, f0, df0, ddf0, dn, a):
    """
    Returns the numerical solutions for the Falkner-Skan equations using
    the fourth-order Runge-Kutta method.
    """
    N = int((nf - ni) / dn)
    n = np.linspace(ni, nf,N+1)
    scaled = ( ( 0.5 * (a+1) )**(1/2) ) * n
    F = np.zeros(len(n))
    dF = np.zeros(len(n))
    ddF = np.zeros(len(n))
    F[0] = f0
    dF[0] = df0
    ddF[0] = ddf0
    for i in range(len(n)-1):
        k01 = g1(F[i], dF[i], ddF[i], a, n[i])
        k02 = g2(F[i], dF[i], ddF[i], a, n[i])
        k03 = g3(F[i], dF[i], ddF[i], a, n[i])
        
        k11 = g1(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, a, n[i]+dn/2)
        k12 = g2(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, a, n[i]+dn/2)
        k13 = g3(F[i]+dn/2*k01, dF[i]+dn/2*k02, ddF[i]+dn/2*k03, a, n[i]+dn/2)
        
        k21 = g1(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, a, n[i]+dn/2)
        k22 = g2(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, a, n[i]+dn/2)
        k23 = g3(F[i]+dn/2*k11, dF[i]+dn/2*k12, ddF[i]+dn/2*k13, a, n[i]+dn/2)
        
        k31 = g1(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, a, n[i]+dn)
        k32 = g2(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, a, n[i]+dn)
        k33 = g3(F[i]+dn*k21, dF[i]+dn*k22, ddF[i]+dn*k23, a, n[i]+dn)
        
        F[i+1] = F[i] + dn/6 * (k01 + 2*k11 + 2*k21 + k31)
        dF[i+1] = dF[i] + dn/6 * (k02 + 2*k12 + 2*k22 + k32)
        ddF[i+1] = ddF[i] + dn/6 * (k03 + 2*k13 + 2*k23 + k33)
    return n, F, dF, ddF, scaled

# Define problem-specific variables
f0 = 0         # initial value of the dependent variable f
df0 = 0        # initial value of the first derivative of f
ddf0 = 0.3321  # initial value of the second derivative of f
ni = 0.0       # initial value of the independent variable n
nf = 5.0       # final value of the independent variable n
dn = 0.01      # step size
a  = 0         # constant parameter of the Falkner-Skan equations

# Solve the Falkner-Skan equations using the Runge-Kutta method
n, F, dF, ddF, scaled = rk4(ni, nf, f0, df0, ddf0, dn, a)

# Plot the solution for f'(eta)
ax.plot(dF, n, label='Newtonian')

ax.legend(loc='upper left')
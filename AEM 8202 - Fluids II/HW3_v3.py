# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 21:44:28 2023

@author: jjser
Justine John A. Serdoncillo 
AEM 8202 Fluids II
Homework 3 v3
    - No more functions
"""
import numpy as np
import matplotlib.pyplot as plt

# %%
"""
    Personal Created Functions
    - g1, g2 and g3 represents 3 linear ODEs
    - rk4 is the Runge-Kutta 4th Order Method
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
    return n, F, dF, ddF

# %% 
"""
    Blasius Boundary Condition
    alpha = 0
    shear stress or f''(0) = 0.3321
"""
def blasius():
    # Define problem-specific variables
    f0 = 0        # initial value of the dependent variable f
    df0 = 0       # initial value of the first derivative of f
    ddf0 = 0.3321  # initial value of the second derivative of f
    ni = 0.0      # initial value of the independent variable n
    nf = 10.0      # final value of the independent variable n
    dn = 0.2      # step size
    a = 0         # constant parameter of the Falkner-Skan equations
    
    # Solve the Falkner-Skan equations using the Runge-Kutta method
    n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn, a)
    
    # Plot the numerical solutions for f, f', and f'' as functions of n
    #plt.plot(n, F, label='f(n)')
    plt.plot(n, dF, label="f'(n)")
    plt.plot(n, ddF, label="f''(n)")
    plt.xlabel("\\eta")
    plt.legend()
    plt.show()

# %% 
"""
    Incipient Separation
    alpha = -0.0904
    shear stress or f''(0) = 0
"""
def sep():
    # Define problem-specific variables
    f0 = 0        # initial value of the dependent variable f
    df0 = 0       # initial value of the first derivative of f
    ddf0 = 0.0  # initial value of the second derivative of f
    ni = 0.0      # initial value of the independent variable n
    nf = 10.0      # final value of the independent variable n
    dn = 0.2      # step size
    a = -0.0904         # constant parameter of the Falkner-Skan equations
    
    # Solve the Falkner-Skan equations using the Runge-Kutta method
    n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn, a)
    
    # Plot the numerical solutions for f, f', and f'' as functions of n
    #plt.plot(n, F, label='f(n)')
    plt.plot(n, dF, label="f'(n)")
    #plt.plot(n, ddF, label="f''(n)")
    plt.xlabel("$ \\eta $")
    plt.legend()
    plt.show()
    
# %% 
"""
    Acceleration
    alpha = 1/3
    shear stress or f''(0) = 0.7570
"""
# Define problem-specific variables
f0 = 0        # initial value of the dependent variable f
df0 = 0       # initial value of the first derivative of f
ddf0 = 0.757  # initial value of the second derivative of f
ni = 0.0      # initial value of the independent variable n
nf = 10.0      # final value of the independent variable n
dn = 0.001      # step size
a = 1/3         # constant parameter of the Falkner-Skan equations

# Solve the Falkner-Skan equations using the Runge-Kutta method
n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn, a)

# Plot the numerical solutions for f, f', and f'' as functions of n
#plt.plot(n, F, label='f(n)')
plt.plot(n, dF, label="f'(n)")
plt.plot(np.linspace(0,10,100,True), np.ones(100), '--')
#plt.plot(n, ddF, label="f''(n)")
plt.xlabel("$ \\eta $")
plt.legend()
plt.show()

# %% 
"""
    Gradual Acceleration
    alpha = 1/9
    shear stress or f''(0) = 0.512
"""
def grad():
    # Define problem-specific variables
    f0 = 0        # initial value of the dependent variable f
    df0 = 0       # initial value of the first derivative of f
    ddf0 = 0.512  # initial value of the second derivative of f
    ni = 0.0      # initial value of the independent variable n
    nf = 10.0      # final value of the independent variable n
    dn = 0.001      # step size
    a = 1/9         # constant parameter of the Falkner-Skan equations
    
    # Solve the Falkner-Skan equations using the Runge-Kutta method
    n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn, a)
    
    # Plot the numerical solutions for f, f', and f'' as functions of n
    #plt.plot(n, F, label='f(n)')
    plt.plot(n, dF, label="f'(n)")
    plt.plot(np.linspace(0,10,100,True), np.ones(100), '--')
    #plt.plot(n, ddF, label="f''(n)")
    plt.xlabel("$ \\eta $")
    plt.legend()
    plt.show()

# %% 
"""
    Linear Acceleration
    alpha = 1
    shear stress or f''(0) = 1.2325
"""
def lin():
    # Define problem-specific variables
    f0 = 0        # initial value of the dependent variable f
    df0 = 0       # initial value of the first derivative of f
    ddf0 = 1.2325  # initial value of the second derivative of f
    ni = 0.0      # initial value of the independent variable n
    nf = 10.0      # final value of the independent variable n
    dn = 0.001      # step size
    a = 1         # constant parameter of the Falkner-Skan equations
    
    # Solve the Falkner-Skan equations using the Runge-Kutta method
    n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn, a)
    
    # Plot the numerical solutions for f, f', and f'' as functions of n
    #plt.plot(n, F, label='f(n)')
    plt.plot(n, dF, label="f'(n)")
    plt.plot(np.linspace(0,10,100,True), np.ones(100), '--')
    #plt.plot(n, ddF, label="f''(n)")
    plt.xlabel("$ \\eta $")
    plt.legend()
    plt.show()
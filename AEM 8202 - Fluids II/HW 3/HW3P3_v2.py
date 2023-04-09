# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 21:44:28 2023

@author: jjser
Justine John A. Serdoncillo 
AEM 8202 Fluids II
Homework 3 Problem 3 v2
    - Solve the Non-Newtonian Similarity ODE
    - Working shooting method
    - added shooting
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings

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
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = -f * ddf / (ddf**(1/2))
            if w:
                print(f"Warning occurred with: \n f={np.around(f,6)},\n df={np.around(df,6)},\n ddf={np.around(ddf,6)},\n n={np.around(n,6)}")
                print(f"Warning message: {w[0].message}")
    except ZeroDivisionError:
        print(f"ZeroDivisionError occurred with  \n f={np.around(f,6)},\n df={np.around(df,6)},\n ddf={np.around(ddf,6)},\n n={np.around(n,6)}")
        result = 0  # or set a suitable value when a division by zero occurs
    return result

def rk4(ni, nf, f0, df0, ddf0, dn):
    N = int((nf - ni) / dn)
    n = np.linspace(ni, nf, N+1)
    F = np.ones(len(n))
    dF = np.ones(len(n))
    ddF = np.ones(len(n))
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
        
        if dF[i+1] > 3:
            #print("bruh you high messed up")
            dF[i+1] = 3
            break
        elif dF[i+1] < -1:
            #print("bruh you low messed up")
            dF[i+1] = 0
            break
        
    return n, F, dF, ddF

def shoot(ni, nf, f0, df0, dn, up, down, show = False):
    ite = 0 
    MAXITE = 100
    bad = True
    TOL = 10E-14
    target = 1
    while bad and ite < MAXITE:
        skrt = (up + down)/2
        ite += 1
        n, F, dF, ddF = rk4(ni, nf, f0, df0, skrt, dn)
        if np.max(dF) > target + TOL:
            down = skrt
        elif np.max(dF) < target - TOL:
            up = skrt
        else:
            bad = False
        if show:
            print(f"I shot with {skrt}")
            print(f"It's been {ite} days since cookies")
            print(f"This is my {np.around(np.max(dF),6)} power \n ")
            #print(f"For the last time {np.around(dF[-1],6)} ")
            plt.plot(n, dF, label=str(ite))
    dff0 = skrt
    if show:
        plt.title(f"This is the {ite} time bro")
        plt.legend()
    return dff0, ite

#yay, skrt = shoot(0, 2.5, 0, 0, 0.001, 0.5, 0.75, True)
## using (0.5,0.75) and 10E-14, (-1, 3)
# -> 0.618964749629 in 41 days

# %% 
"""
   Non-Newtonian Boundary Layer Profile
   - dddf * ddf^(1/2) + f * ddf = 0
"""

# Define problem-specific variables
f0 = 0        # initial value of the dependent variable f
df0 = 0       # initial value of the first derivative of f
ni = 0.0      # initial value of the independent variable n
nf = 2.5      # final value of the independent variable n
dn = 0.001     # step size

# Calculate wall shear stress using shooting method
#ddf0, ite = shoot(ni, nf, f0, df0, dn, 0.5, 0.75)
ddf0 = 0.619  # initial value of the second derivative of f

# Solve the Falkner-Skan equations using the Runge-Kutta method
n, F, dF, ddF = rk4(ni, nf, f0, df0, ddf0, dn)
print(f"For Non-Newtonian Boundary Layers, ddf(0) = {np.around(ddf0,4)}")

# Create a plot to display the solution
fig, ax = plt.subplots(figsize=(7,5))

# Set the title and axis labels for the plot
ax.set_title('Comparison betweeen Newtonian and Non-Newtonian \n Boundary Layers')
ax.set_ylabel('$ \\eta $')
ax.set_xlabel('$ f\'(\\eta) $')

# Set the x and y limits for the plot
ax.set_ylim([0,7.0])
ax.set_xlim([0,1.0])

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
ni = 0.0       # initial value of the independent variable n
nf = 10.0      # final value of the independent variable n
dn = 0.001     # step size
a  = 0         # constant parameter of the Falkner-Skan equations

ddf0 = 0.3321  # initial value of the second derivative of f

# Solve the Falkner-Skan equations using the Runge-Kutta method
n, F, dF, ddF, scaled = rk4(ni, nf, f0, df0, ddf0, dn, a)
print(f"For Newtonian Boundary Layers, ddf(0) = {np.around(ddf0,4)}")

# Plot the solution for f'(eta)
ax.plot(dF, n, label='Newtonian')

ax.legend(loc='upper left')

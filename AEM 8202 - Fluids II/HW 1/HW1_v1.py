# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 10:06:00 2023

@author: jjser
"""
# Imported from CFD Project
import numpy as np
import matplotlib.pyplot as plt
import math
from numpy import linalg as la
from numpy import asarray
import pandas
import time

##################
def G(r):
    return np.log(r)/(2*np.pi)
def dGdn(R,N,r):
    return R @ N / (2*np.pi*r)

def prob(a,b,n):
    if a == b:
        shape = 'Circle Cylinder'
    else:
        shape = 'Ellipse Cylinder'
    dtheta = 2*np.pi / n
    xO = np.zeros(n+1)
    yO = np.zeros(n+1)
    theta = np.zeros(n+1)
    
    # Creating the actual circle first
    for i in range(n+1):
        theta[i] = i * dtheta
        r = math.sqrt((a*b)/(a**2*np.sin(theta[i])**2+b**2*np.cos(theta[i])**2))
        xO[i] = r*np.cos(theta[i])
        yO[i] = r*np.sin(theta[i])
    # Applying BEM
    A = np.zeros((n,n))
    B = np.zeros((n,n))
    area = np.zeros(n)
    cx = np.zeros(n)
    cy = np.zeros(n)
    N = np.zeros((n,2))
    
    s = 2
    fig0, ax0 = plt.subplots()
    title = 'BEM for ' + shape + ', n = ' + str(n)
    ax0.set_title(title)
    ax0.set_xlabel('x')
    ax0.set_ylabel('y')
    ax0.set_xlim([-s,s])
    ax0.set_ylim([-s,s])
    
    for i in range(n):
        area[i] = math.sqrt((xO[i+1]-xO[i])**2+(yO[i+1]-yO[i])**2)
        cx[i] = 0.5 *  (xO[i+1]+xO[i])
        cy[i] = 0.5 *  (yO[i+1]+yO[i])
        N[i][0] = -(yO[i+1]-yO[i])
        N[i][1] =  (xO[i+1]-xO[i])
        N[i] = N[i] / la.norm(N[i])
    for i in range(n):
        for j in range(n):
            if i==j:
                B[i,j] = 0
                A[i,j] = -0.5
            else:
                R = [cx[j]-cx[i],cy[j]-cy[i]]
                r = la.norm(R)
                R = R/r
                A[i,j] = dGdn(R,N[j],r) * area[j]
                B[i,j] = G(r) * area[j]
                #ax0.quiver(cx[i],cy[i],R[0],R[1])
                
    # plotting 
    ax0.plot(xO,yO)
    #ax0.plot(cx,cy,'r*')
    #ax0.quiver(cx,cy,N[:,0],N[:,1])
    #ax0.contour(x, y, psi)
    # psi1 and dpsi1dn
    psi1 = np.zeros(n)
    for i in range(n):
        psi1[i] = -yO[i]
    #print(psi1)
    #print(area)
    #print(B)
    #print(A)
    
    dpsi1dn = la.inv(B) @ A @ psi1.T
    
    # Plot psi1 dpsi1dn
    fig, ax = plt.subplots()
    title = '$\\psi_1$ and $\\frac{\\partial \\psi_1}{n_1}$ for ' + shape + ', n = ' + str(n) 
    ax.set_title(title)
    ax.set_xlabel('$\\theta$')
    
    ax.plot(theta[:-1],psi1,label='$\\psi_1$')
    ax.plot(theta[:-1],dpsi1dn,label='$\\frac{\\partial \\psi_1}{n_1}$')
    ax.legend()
    if a > b:
        #ax.set_ylim([-1.25,1.25])
        pass
    else:
        #ax.set_ylim([-2,2])
        pass
    
    # Plotting the Cylinder
    g = 100
    x = np.linspace(-s,s,g)
    y = np.linspace(-s,s,g)
    X,Y = np.meshgrid(x,y)
    psi = np.zeros((g,g))
    '''
    for i in range(g):
        for j in range(g):
            for k in range(n):
                rx = cx[k] - X[i,j]
                ry = cy[k] - Y[i,j]
                rmag = np.sqrt(rx**2+ry**2)
                R = np.asarray([rx,ry]/rmag) 
                psi += (R@N[k]/(2*np.pi*rmag) - np.log(rmag)/(2*np.pi)*dpsi1dn[k])*area[k]
    psi += Y
    '''
    '''      
    for k in range(n):
        rx = cx[k] - X
        ry = cy[k] - Y 
        rmag = np.sqrt(rx**2+ry**2)
        R = np.asarray([rx,ry]/rmag) 
        #print(R.shape)
        psi += dGdn(R,N[k],rmag) * area[k] - G(rmag) * dpsi1dn[k] * area[k]
    psi += Y
    '''
    

##### BEM for a 4:1 ellipse
prob(1,1,200)
#prob(1,1,100)
#prob(1,1,200)
#prob(1,4,50)
#prob(1,4,100)
#prob(1,4,200)

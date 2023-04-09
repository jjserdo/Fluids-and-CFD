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
from scipy.stats import norm 

##################
def G(r):
    return np.log(r)/(2*np.pi)
def dGdn(R,N,r):
    return R @ N / (2*np.pi*r)

def prob(a,b,n,num,prob):
    # 'a' is the length of the major axis
    # 'b' is the length of the minor axis
    # 'n' is the number of elements in the boundary element method
    # 'num' is the number of streamlines
    # 'prob' is either 
    if a > b:
        shape = 'Greenhouse'
        prob = 1
    elif b > a:
        shape = 'Prized Rose'
        prob = 2
    else:
        shape = 'bro wrong problem'
        prob = False
        
    dtheta = 2*np.pi / n
    xO = np.zeros(n+1)
    yO = np.zeros(n+1)
    theta = np.zeros(n+1)
    
    # Creating the actual circle first
    for i in range(n+1):
        theta[i] = i * dtheta
        xO[i] = a*np.cos(theta[i])
        yO[i] = b*np.sin(theta[i])
    # Applying BEM
    A = np.zeros((n,n))
    B = np.zeros((n,n))
    area = np.zeros(n)
    cx = np.zeros(n)
    cy = np.zeros(n)
    N = np.zeros((n,2))
    
    s = 1.5*max(a,b)
    #fig0, ax0 = plt.subplots(figsize=(5,5))
    fig0, ax0 = plt.subplots()
    title = 'BEM for ' + shape + ', n = ' + str(n)
    ax0.set_title(title)
    ax0.set_xlabel('x')
    ax0.set_ylabel('y')
    ax0.set_xlim([-s,s])
    ax0.set_ylim([0,1.1*b])
    
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
    ax0.plot(xO,yO)
    
    # Solve for psi1 and dpsi1dn
    psi1 = np.zeros(n)
    for i in range(n):
        psi1[i] = -yO[i]
    dpsi1dn = la.inv(B) @ A @ psi1.T
    
    # Plotting the Cylinder
    g = 501
    x = np.linspace(-s,s,g)
    y = np.linspace(-s,s,g)
    X,Y = np.meshgrid(x,y)
    psi = np.zeros((g,g))
    
    # Streamfunction Contours       
    for k in range(n):
        rx = cx[k] - X
        ry = cy[k] - Y 
        rmag = np.sqrt(rx**2+ry**2)
        psi += (psi1[k]*(rx*N[k][0]+ry*N[k][1])/(2*np.pi*rmag)-G(rmag)*dpsi1dn[k])*area[k]
    psi += Y

    mean = np.mean(psi)
    std = np.std(psi)
    levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
    ax0.contour(x, y, psi, levels=levels)
    
    if prob != False:
        name = 'images/BEM' + str(prob) + '.png'
        fig0.savefig(name)
    
##### BEM for a 4:1 ellipse
#prob(4,1,500,15)
prob(96,48,500,15)
prob(48,96,500,15)
start_time = time.time()
print("--- %10s seconds ---" % np.round((time.time() - start_time),4))

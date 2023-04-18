# -*- coding: utf-8 -*-
"""
Created using prof's notes in Tuesday 12/6 class

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt

def init():
    # Initializing Everything
    ni = 5 # number of cells in i
    nj = 5 # number of cells in j
    g = 1 # number of ghost cells
    x = np.zeros((ni+2*g,nj+2*g))
    y = np.zeros((ni+2*g,nj+2*g))
    Jdet = np.zeros((ni+2*g,nj+2*g))
    dxidx = np.zeros((ni+2*g,nj+2*g))
    dxidy = np.zeros((ni+2*g,nj+2*g))
    detadx = np.zeros((ni+2*g,nj+2*g))
    detady = np.zeros((ni+2*g,nj+2*g))
    
    # Notes: ghost cells are treated as 0 index
    ib = g
    jb = g
    ie = ib + ni - 1
    je = jb + nj - 1
    return ni,nj,g,x,y,Jdet,dxidx,dxidy,detadx,detady,ib,jb,ie,je

def normgrid(ni,nj,x,y,ib,jb,ie,je):
    # Initializing the Grid
    for i in np.arange(ib,ie+1): # goes from ib to ie
        for j in np.arange(jb,je+1): # goes from jb to je
            x[i,j] = (i-1) * 2.0 / (ni-1)
            y[i,j] = (j-1) * 1.0 / (nj-1)
    return x, y

def calc(ni,nj,x,y,Jdet,dxidx,dxidy,detadx,detady,ib,jb,ie,je):
    for i in np.arange(ib,ie):
        for j in np.arange(jb,je):
            if i >=2 and i<= ni-1:
                # central
                dx_dxi = (x[i+1,j]-x[i-1,j]) / 2.0 # 2.0 is 2*dxi
                dy_dxi = (x[i+1,j]-x[i-1,j]) / 2.0
            elif i == 1:
                # forward
                dx_dxi = (x[i+1,j]-x[i,j]) / 1.0
                dy_dxi = (x[i+1,j]-x[i,j]) / 1.0
            elif i == ni:
                # forward
                dx_dxi = (x[i,j]-x[i-1,j]) / 1.0
                dy_dxi = (x[i,j]-x[i-1,j]) / 1.0
            
            if j >=2 and j<= nj-1:
                # central
                dx_deta = (x[i,j+1]-x[i,j-1]) / 2.0 # 2.0 is 2*dxi
                dy_deta = (x[i,j+1]-x[i,j-1]) / 2.0
            elif j == 1:
                # forward
                dx_deta = (x[i,j+1]-x[i,j]) / 1.0
                dy_deta = (x[i,j+1]-x[i,j]) / 1.0
            elif j == nj:
                # forward
                dx_deta = (x[i,j]-x[i,j-1]) / 1.0
                dy_deta = (x[i,j]-x[i,j-1]) / 1.0    
            Jdet[i,j] = 1/(dx_dxi*dy_deta-dx_deta*dy_dxi)
            # note sure
            dxidx[i,j] = dy_deta * Jdet[i,j]
            dxidy[i,j] = -dx_deta * Jdet[i,j]
            detadx[i,j] = -dy_dxi * Jdet[i,j]
            detady[i,j] = dx_dxi * Jdet[i,j]
            
        return Jdet,dxidx,dxidy,detadx,detady

def main():
    ni,nj,g,x,y,Jdet,dxidx,dxidy,detadx,detady,ib,jb,ie,je = init()
    x,y = normgrid(ni,nj,x,y,ib,jb,ie,je)  
    Jdet,dxidx,dxidy,detadx,detady = calc(ni,nj,x,y,Jdet,dxidx,dxidy,detadx,detady,ib,jb,ie,je)    

if __name__ == '__main__':
    main()
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 22:41:15 2023

@author: jjser
Justine John A. Serdoncillo
Potential Flow: Airfoil
Initially created for:
    - Spring 2023 
    - AEM 8202 Fluids II 
    - Final Exam
Updates:
    - 
"""

import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.stats import norm 

# %%
class airfoil():
    def __init__(self, c, tmax, N = 100, alpha = 0):
        """
        Parameters
        ----------
        c : float
            - chord length
            - default is 1
        tmax : float
            - % maximum thickness of the airfoil w/ respect to l
        alpha : float 
            - angle of attack of the flow (deg)
            - default is 0
        """
        # Here I assign the input to the class
        self.c = c
        self.tmax = tmax
        self.alpha = alpha
        self.N = N
        # Here I create the grid
        x = np.linspace(-5, 5, 1000)
        y = np.linspace(-5, 5, 1000)
        self.X,self.Y = np.meshgrid(x,y)
        self.z = self.X + 1j * self.Y
        # Here I have my plot stuff
        self.plottedZeta = False
        self.plottedZ = False
        self.rotated = False
        self.kuttaed = False
        # Here I run my functions
        self.getM()
        self.getPoints()
        self.compute()
        self.thick()
        
    def compute(self):
        # Solving for FZeta and FZ
        if self.rotated == True:
            self.FZeta = (self.z * np.exp(-1j*self.alpha) + self.m) + (1+self.m)**2/(self.z * np.exp(-1j*self.alpha) + self.m)
            if self.kuttaed == True:
                self.FZeta += 1j * (1 + self.m) * np.log(self.z)
        else:
            self.FZeta = (self.z + self.m) + (1+self.m)**2/(self.z + self.m)
        self.psiZeta = self.FZeta.imag
        
        self.l = 1/2 * (self.z + np.sqrt(self.z+2) * np.sqrt(self.z-2) )
        if self.rotated == True:
            self.FZ = (self.l * np.exp(-1j*self.alpha) + self.m) + (1+self.m)**2/(self.l * np.exp(-1j*self.alpha) +self.m)
            if self.kuttaed == True:
                self.FZ += 1j * (1 + self.m) * np.log(self.z)
        else:
            self.FZ = (self.l+self.m) + (1+self.m)**2/(self.l+self.m)
        self.psiZ = self.FZ.imag
        
    def plotZeta(self):
        if not self.plottedZeta:
            self.figZeta, self.axZeta = plt.subplots(figsize=(10,10)) 
            self.plottedZeta = True
        self.axZeta.set_title("Plotted Zeta Contour")
        self.axZeta.set_xlim([-1.5, 1.5])
        self.axZeta.set_ylim([-1.5, 1.5])
        mean = np.mean(self.psiZeta)
        std = np.std(self.psiZeta)
        levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), 20) 
        self.axZeta.contour(self.X, self.Y, self.psiZeta, levels=levels)    
        
        return True
    
    def plotZ(self):
        if not self.plottedZ:
            self.figZ, self.axZ = plt.subplots(figsize=(10,10))
            self.plottedZ = True    
        self.axZ.set_title("Plotted Z Contour")
        self.axZ.set_xlim([-2.5, 2.5])
        self.axZ.set_ylim([-1.5, 1.5])
        #self.axZ.set_ylim([-1.5 * self.thickness, 1.5 * self.thickness])
        mean = np.mean(self.psiZ)
        std = np.std(self.psiZ)
        levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), 100) 
        self.axZ.contour(self.X, self.Y, self.psiZ, levels=levels)  

        return True
    
    def plotPointsZeta(self):
        if not self.plottedZeta:
            figZeta, axZeta = plt.subplots(figsize=(10,10))
            self.plottedZeta = True
        
        self.axZeta.set_title("Plotted Zeta Points")
        self.axZeta.plot(self.OGZeta.real, self.OGZeta.imag, 'C0')
        self.axZeta.plot(self.OGZeta.real, -self.OGZeta.imag, 'C0')
        self.axZeta.plot(self.pointsZeta.real, self.pointsZeta.imag, 'C1')
        self.axZeta.plot(self.pointsZeta.real, -self.pointsZeta.imag, 'C1')
        
        return True
    
    def plotPointsZ(self):
        if not self.plottedZ:
            self.figZ, self.axZ = plt.subplots(figsize=(10,10))
            self.plottedZ = True
        
        self.axZ.set_title("Plotted Z Points")
        self.axZ.plot(self.OGZ.real, self.OGZ.imag, 'C0')
        self.axZ.plot(self.pointsZ.real, self.pointsZ.imag, 'C1')
        self.axZ.plot(self.pointsZ.real, -self.pointsZ.imag, 'C1')
        
        return True
        
    def trans(self, zetaSpace):
        zSpace = zetaSpace + 1/zetaSpace
        return zSpace
    
    def getChordLength(self):
        self.chordLength = np.max(self.pointsZ.real) - np.min(self.pointsZ.real)
        return True
    
    def getPoints(self):
        # Point Properties
        self.thetas = np.linspace(np.pi, 0, self.N+1, endpoint=True)
        self.OGZeta = np.cos(self.thetas) + 1j * np.sin(self.thetas)
        self.pointsZeta = ( (1+self.m)*np.cos(self.thetas) - self.m ) + 1j * (1+self.m) * np.sin(self.thetas)
        self.OGZ = self.trans(self.OGZeta)
        self.pointsZ = self.trans(self.pointsZeta)
        self.dsCum = np.zeros(self.N+1)

        # Line Properties
        self.ds = np.zeros((self.N,2))
        self.sXmid = np.zeros(self.N)
        self.sYmid = np.zeros(self.N)
        self.dsZ = np.zeros((self.N,2))
        self.sXmidZ = np.zeros(self.N)
        self.sYmidZ = np.zeros(self.N)
        
        
        for ii in range(self.N):
            self.ds[ii] = [self.pointsZeta.real[ii+1]-self.pointsZeta.real[ii],self.pointsZeta.imag[ii+1]-self.pointsZeta.imag[ii]]
            self.sXmid[ii] = 0.5 * (self.pointsZeta.real[ii+1]+self.pointsZeta.real[ii])
            self.sYmid[ii] = 0.5 * (self.pointsZeta.imag[ii+1]+self.pointsZeta.imag[ii])
            self.dsZ[ii] = [self.pointsZ.real[ii+1]-self.pointsZ.real[ii],self.pointsZ.imag[ii+1]-self.pointsZ.imag[ii]]
            self.sXmidZ[ii] = 0.5 * (self.pointsZ.real[ii+1]+self.pointsZ.real[ii])
            self.sYmidZ[ii] = 0.5 * (self.pointsZ.imag[ii+1]+self.pointsZ.imag[ii])
            
        self.dsNorm = la.norm(self.ds, axis=1, keepdims=True)
        self.dsDir = self.ds/self.dsNorm
        self.dsNormZ = la.norm(self.dsZ, axis=1, keepdims=True)
        self.dsDirZ = self.ds/self.dsNormZ
        
        for ii in range(self.N):
            self.dsCum[ii+1] = self.dsCum[ii]
            self.dsCum[ii+1] += self.dsNormZ[ii]
            
    def thick(self):
        maxThicc = np.max(self.pointsZ.imag)
        self.thickness = 2 * maxThicc
        
    def getM(self):
        flag = True # set to True for Speed
        TOL = 10E-8
        ite = 0
        MAXITE = 100
        up = 0.5
        down = 0
        diff = 0
        
        while not flag and ite < MAXITE:
            self.m = (up + down)/2
            ite += 1
            self.getPoints()
            self.compute()
            self.plotZ()
            maxThicc = np.max(self.pointsZ.imag)
            #print(f"I tried {ite} times")
            #print(f"I am m right now {self.m}")
            #print(f"I am thickness {2*maxThicc}")
            self.chordLength = np.max(self.pointsZ.real) - np.min(self.pointsZ.real)
            #print(f"I am chord length {self.chordLength}")
            #print(f"I am LEGEND {self.chordLength * self.tmax/100}")
            diff = maxThicc * 2 - self.chordLength * self.tmax/100
            #print(f"we differ by {diff}")
            if np.abs(diff) > TOL:
                if diff > 0:
                    up = self.m
                else:
                    down = self.m
            else:
                self.thickness = 2 * maxThicc
                flag = True
        #print(ite)
        #print(diff)
        if self.tmax == 8:
            self.m = 0.0656920075416565
        elif self.tmax == 4:
            self.m = 0.03177905082702637
        elif self.tmax == 2:
            self.m = 0.01563858985900879
        
        return True
    
    def wwStar(self, zeta, z):
        w = (1 - (1 + self.m)**2/(zeta + self.m)**2) * (1/2 * (1 + z / (np.sqrt(z+2) * np.sqrt(z-2))))
        wStar = np.conj(w)
        return w, wStar
    
    def plotTangentialVelocities(self):
        # Point Properties
        self.velocities = np.zeros(self.N+1)
        self.wVel = np.zeros(self.N+1)
        for ii in range(self.N+1):
            self.wVel[ii], wStar = self.wwStar(self.pointsZeta[ii], self.pointsZ[ii])
            self.velocities[ii] = self.wVel[ii] * wStar
        #self.speeds = la.norm(self.velocities, axis=1, keepdims=True)
        self.uVel =  self.wVel.real
        self.vVel = -self.wVel.imag
        
        # Line Properties
        self.dUds = np.zeros(self.N)
        for ii in range(self.N):
            self.dUds[ii] = self.velocities[ii+1] - self.velocities[ii]    
        fig, ax = plt.subplots()
        ax.set_title("Tangential Velocity on the top surface")
        ax.plot(self.dsCum, self.velocities)
        
        return True
    
    def rotateFlow(self, alpha = 0, kutta = False):
        """
        Parameters
        ----------
        alpha : float
            - new angle of attack (deg)
        """
        # Calculate the Circulation required for Kutta Condition
        if alpha != 0:
            self.rotated = True
            self.alpha = np.pi * alpha/180
            if kutta == True:
                self.kuttaed = True
            self.compute()
            # Plot New Streamlines and Points
            self.plotZeta()
            self.plotPointsZeta()
            self.plotZ()
            self.plotPointsZ()
        else:
            print('Do not fool me')
        
        return True
    
    def getLift_Blasius(self, rho):
        """
        Parameters
        ----------
        rho : float
            - flow density 
        """
        R = 2.5 
        f = 0
        for ii in range(self.N-1):
            ''' Not sure what to put here yet '''
            f = 0     
        f = 1j * rho / 2 * f
        fx = f.real
        fy = f.imag
        self.liftBlasius = fy * np.cos(self.alpha) + fx * np.sin(self.alpha)
        return self.liftBlasius
    
    def doThwaites(self, nu):
        """
        Parameters
        ----------
        nu : float
            - kinematic viscosity
        """
        self.theta = np.zeros(self.N+1)
        self.delta = np.zeros(self.N+1)
        integral = 0
        for ii in range(self.N):
            integral += self.velocities[ii]**5 * self.dsNorm[ii]
            self.theta[ii] = np.sqrt(0.45*nu*integral)
            self.delta[ii] = self.theta**2/nu*self.dUds
        fig, ax = plt.subplots()
        ax.set_title('Thwaites $ \\theta $')
        ax.plot(self.dsCum, self.theta)
        
        fig1, ax1 = plt.subplots()
        ax.set_title('Thwaites $ \\delta $')
        ax1.plot(self.dsCum, self.delta)
        return True
    

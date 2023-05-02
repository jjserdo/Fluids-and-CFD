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
import matplotlib.pyplot as plt

class airfoil():
    def __init__(self, c, tmax, alpha = 0):
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
        self.c = c
        self.tmax = tmax
        self.alpha = alpha
        self.getM()
        
    def getM(self):
        self.m = self.tmax
        return True
    
    def compute(self):
        
        return True
    
    def plotStreamlines(self):
        
        return True
    
    def plotTangentialVelocities(self):
        
        return True
    
    def rotateFlow(self, alpha):
        """
        Parameters
        ----------
        alpha : float
            - new angle of attack (deg)
        """
        self.compute()
        return True
    
    def getLift_Blasius(self, rho):
        """
        Parameters
        ----------
        rho : float
            - flow density 
        """
        lift = 0
        return lift
    
    def getLift_Integrate(self, rho):
        """
        Parameters
        ----------
        rho : float
            - flow density 
        """
        lift = 0
        return lift
    
    def doThwaites(self, nu):
        """
        Parameters
        ----------
        nu : float
            - kinematic viscosity
        """
        theta = 0
        delta = 0
        return theta, delta
    

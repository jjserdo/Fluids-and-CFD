# Imported from CFD Project
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import norm 
##################
def HW4(problem='blah',let='a',num=None,num2=None,
        plotzero=False,plotstag=False,func='nope'):
    s = 5
    g = 500
    X = np.linspace(-s,s,g)
    Y = np.linspace(-s,s,g)
    
    fig, ax = plt.subplots(figsize=(5,5))
    #ax.grid()
    ax.set_title('Complex Plots')
    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    ax.set_xlim([-s,s])
    ax.set_ylim([-s,s])
    
    if problem == 'Problem 1':
        pass    
    elif problem == 'Problem 2':
        s = 5
        X = np.linspace(-s,   s, g)
        Y = np.linspace(-s, 2*s, g)
        ax.set_xlim([-2,  2])
        ax.set_ylim([ 0,  2])
    elif problem == 'Problem 3':
        if let == 'b':
            s = 2
            X = np.linspace(-192,192,g)
            Y = np.linspace(   0, 70,g)
            ax.set_xlim([-192,192])
            ax.set_ylim([   0, 60])
        if let == 'c':
            s = 2
            X = np.linspace(-96, 96,g)
            Y = np.linspace(  0,130,g)
            ax.set_xlim([-96, 96])
            ax.set_ylim([  0,130])
        
    x,y = np.meshgrid(X,Y)
    r = np.arctan2(y,x)
    t = np.sqrt(x**2+y**2)
    z = x + 1j*y
    ''' Input Function Here '''
    U = 2/3
    a = np.pi/3
    c = U*np.exp(-1j*a)
    Q = 5
    G = 1
    M = 1
    A = 1
    sep = 2
    
    func = 'lmao'
    if problem == 'Problem 1':
        Fz = -1j*G/(2*np.pi)*np.log(z**3-1) + U*z**(1.5)
    elif problem == 'Problem 2':
        Fz =  1j*M/(c*(z-1j)) + z**2
    elif problem == 'Problem 3':
        if let == 'b': 
            U = 1
            Q = 1
            sep = 1
        if let == 'c':
            U = 1
            Q = 1
            sep = 1
        Fz = U*z + Q/(2*np.pi)*np.log((z+sep)/(z-sep))
        #Fz = Q*sep/(np.pi*z)
    else:
        #Fz = U*(z+A**2/z) + 1j*G/(2*np.pi)*np.log(z)
        # Success
        # z, z**2, c*z, Q/(2*np.pi)*np.log(z), -1j*G/(2*np.pi)*np.log(z)
        # U*z**(1.5), M/z, U*z + M/z, U*(z+A**2/z)
        #Fz = M/z
        Fz = z**4/4-z
        
    phi = Fz.real
    psi = Fz.imag
    u = np.zeros((g,g)) # dpsidy = dphidx
    v = np.zeros((g,g)) # -dpsidx = dphidy
    dx = X[1] - X[0]
    dy = Y[1] - Y[0]
    
    ''' u and v functions are weird? '''
    for i in range(g):
        for j in range(g):
            if i == 0:
                v[i,j] = (phi[i+1,j]-phi[i,j]) / dx
            elif i == g-1:
                v[i,j] = (phi[i,j]-phi[i-1,j]) / dx
            else:
                v[i,j] = (phi[i+1,j]-phi[i-1,j]) / (2*dx)
            if i == 0:
                u[i,j] = (psi[i+1,j]-psi[i,j]) / dx
            elif i == g-1:
                u[i,j] = (psi[i,j]-psi[i-1,j]) / dx
            else:
                u[i,j] = (psi[i+1,j]-psi[i-1,j]) / (2*dx)
    if num != None:
        if func == 'lmao':
            mean = np.mean(phi)
            std = np.std(phi)
            levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
            ax.contour(x, y, phi, levels=levels, linestyles='dashed',alpha=0.5)
            mean = np.mean(psi)
            std = np.std(psi)
            levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
            ax.contour(x, y, psi, levels=levels)
        else:
            ax.contour(x, y, phi, num, linestyles='dashed',alpha=0.5)
            ax.contour(x, y, psi, num)
        
    #### Plot Stagnation Points ####
    mag = np.sqrt(u**2 + v**2)
    if plotstag == True:
        zero = np.argwhere(np.around(mag,1) == 0)
        zeros = np.zeros((2,zero.shape[0]))
        for i in range(zero.shape[0]):
            zeros[0,i] = x[zero[i,0],zero[i,1]]
            zeros[1,i] = y[zero[i,0],zero[i,1]]
        ax.plot(zeros[0],zeros[1],'r*')
        
    #### Plot Zero Streamline ####
    if plotzero == True:
        zero = np.argwhere(np.around(psi,1) == 0)
        zeros = np.zeros((2,zero.shape[0]))
        for i in range(zero.shape[0]):
            zeros[0,i] = x[zero[i,0],zero[i,1]]
            zeros[1,i] = y[zero[i,0],zero[i,1]]
        ax.plot(zeros[0],zeros[1],'b*')
    if num2 != None:
        ll = int(g/num2)
        skip = (slice(None, None, ll), slice(None, None, ll))
        ax.quiver(x[skip],y[skip],u[skip],v[skip])
if __name__ == "__main__":
    #W4('Problem 1',num=10,plotstag=1)
    HW4('Problem 2',func='lmao')
    #HW4('Problem 3','b')
    #HW4('Problem 3','c')
    #HW4(num=10,plotstag=True)


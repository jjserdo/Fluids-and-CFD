# Imported from CFD Project
import numpy as np
import time
from scipy.stats import norm 
import matplotlib.pyplot as plt
##################
def HW4(problem='blah',let='a',num=10,num2=None,
        plotzero=False,plotstag=False,func='nope',s=5):
    g = 500
    X = np.linspace(-s,s,g)
    Y = np.linspace(-s,s,g)
    
    #fig, ax = plt.subplots(figsize=(5,5))
    fig, ax = plt.subplots()
    #ax.grid()
    ax.set_title('Complex Plots')
    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    ax.set_xlim([-s,s])
    ax.set_ylim([-s,s])
    
    if problem == 'Problem 1':
        pass    
    elif problem == 'Problem 2':
        ax.set_xlim([-s,s])
        ax.set_ylim([ 0,s])
    elif problem == 'Problem 3':
        if let == 'b':
            X = np.linspace(-192,192,g)
            Y = np.linspace(   0, 70,g)
            ax.set_xlim([-192,192])
            ax.set_ylim([   0, 60])
        if let == 'c':
            X = np.linspace(-96, 96,g)
            Y = np.linspace(  0,130,g)
            ax.set_xlim([-96, 96])
            ax.set_ylim([  0,130])
        
    x,y = np.meshgrid(X,Y)
    r = np.arctan2(y,x)
    t = np.sqrt(x**2+y**2)
    z = x + 1j*y
    ''' Input Function Here '''
    
    if problem == 'Problem 1':
        Fz = z**4/4 - z
        w = z**3 - 1
        u = w.real
        v = w.imag
    elif problem == 'Problem 2':
        M = 1
        Fz = M*(np.sqrt(3)-z) / (z**2+1)
        w = (z**2 -2*np.sqrt(3)*z-1)/ (z**2+1)**2
        u = w.real
        v = w.imag
    elif problem == 'Problem 3':
        if let == 'b': 
            a = 72
            b = 24 * np.sqrt(3)
        if let == 'c':
            a = 72
            b = 24 * np.sqrt(3) * 1j
        Fz = z/2 * (1+a**2/b**2) + (1-a**2/b**2) * np.sqrt((z/2)**2 - b**2)
    else:
        U = 1
        Q = 10
        sep = 1
        Fz = U*z + Q/(2*np.pi)*np.log((z+sep)/(z-sep))
        
    phi = Fz.real
    psi = Fz.imag
                
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
        
    #### Plot Zero Streamline ####
    if plotzero == True:
        zero = np.argwhere(np.around(psi,1) == 0)
        zeros = np.zeros((2,zero.shape[0]))
        for i in range(zero.shape[0]):
            zeros[0,i] = x[zero[i,0],zero[i,1]]
            zeros[1,i] = y[zero[i,0],zero[i,1]]
        ax.plot(zeros[0],zeros[1],'b*')
        
    if 'u' and 'v' in locals():
        mag = np.sqrt(u**2+v**2)
        #### Plot Stagnation Points ####
        if plotstag == True:
            zero = np.argwhere(np.around(mag,1) == 0)
            zeros = np.zeros((2,zero.shape[0]))
            for i in range(zero.shape[0]):
                zeros[0,i] = x[zero[i,0],zero[i,1]]
                zeros[1,i] = y[zero[i,0],zero[i,1]]
            ax.plot(zeros[0],zeros[1],'r*')
        if num2 != None:
            ll = int(g/num2)
            skip = (slice(None, None, ll), slice(None, None, ll))
            ax.quiver(x[skip],y[skip],u[skip],v[skip])

if __name__ == "__main__":
    #HW4('Problem 1',func='lmao',num=10,plotstag=True,num2=20)
    #HW4('Problem 2',func='lmao',num=10,s=10,plotstag=True) # weird stagnations stuff
    HW4('Problem 3','b',func='lmao')
    #HW4('Problem 3','c',func='lmao')
    #HW4(func='lmao')


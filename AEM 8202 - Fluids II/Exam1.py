# Imported from CFD Project
import numpy as np
import time
from scipy.stats import norm 
import matplotlib.pyplot as plt
##################
def HW4(problem='blah',let='a',num=10,num2=None,
        plotzero=False,plotstag=False,func='nope',s=5,pot=True,stream=False):
    g = 500
    X = np.linspace(-s,s,g)
    Y = np.linspace(-s,s,g)
    '''
    if problem == 'Problem 1':
        fig, ax = plt.subplots(figsize=(5,5)) 
    elif problem == 'Problem 2':
        fig, ax = plt.subplots(figsize=(10,5)) 
    elif problem == 'Problem 3':
        if let == 'b': 
            fig, ax = plt.subplots(figsize=(12,3)) 
        if let == 'c':
            fig, ax = plt.subplots(figsize=(12,5)) 
    if problem == 'blah':
        fig, ax = plt.subplots(figsize=(20,20))
    '''
    fig, ax = plt.subplots()
    #ax.grid()
    ax.set_title('Complex Plots')
    ax.set_xlabel('Real')
    ax.set_ylabel('Imaginary')
    ax.set_xlim([-s,s])
    ax.set_ylim([-s,s])
    
    if problem == 'Problem 1':
        ax.set_xlim([-s,s])
        ax.set_ylim([ 0,s])   
    elif problem == 'Problem 2':
        ax.set_xlim([0,s])
        ax.set_ylim([0,s]) 
    elif problem == 'Problem 3':
        ax.set_xlim([-s,s])
        ax.set_ylim([-s,s]) 
        
    x,y = np.meshgrid(X,Y)
    r = np.arctan2(y,x)
    t = np.sqrt(x**2+y**2)
    z = x + 1j*y
    ''' Input Function Here '''
    
    if problem == 'Problem 1':
        Fz = -1j*np.log((z-2*1j)/(z+2*1j))
        w = -1j/(z-2*1j) + 1j/(z+2*1j)
        u =  w.real
        v = -w.imag
    elif problem == 'Problem 2':
        Fz = -1j*np.log((z**2-2*1j)/(z**2+2*1j))
        w = -1j*2*z/(z**2-2*1j) + 1j*2*z/(z**2+2*1j)
        u =  w.real
        v = -w.imag
    elif problem == 'Problem 3':
        #Fz = -1j*np.log((z**1.5-2*1j)/(z**1.5+2*1j))
        Fz = -1j*np.log(((z*np.exp(-1j*np.pi/3))**1.5*np.exp(1j*np.pi/2)-2*1j)/((z*np.exp(-1j*np.pi/3))**1.5*np.exp(1j*np.pi/2)+2*1j))
        #w = -1j*1.5*z/(z**1.5-2*1j) + 1j*1.5*z/(z**1.5+2*1j)
        #u =  w.real
        #v = -w.imag
    else:
        #Fz = -1j*np.log(z)
        Fz = -1j*np.log(z-2*1j) + 1j*np.log(z+2*1j)
        w = -1j/(z-2*1j) + 1j/(z+2*1j)
        u =  w.real
        v = -w.imag
        
        
    phi = Fz.real
    psi = Fz.imag
                
    if func == 'lmao':
        if pot == True:
            mean = np.mean(phi)
            std = np.std(phi)
            levels = np.linspace(norm(mean,std).ppf(0.05), norm(mean,std).ppf(0.95), num)  
            ax.contour(x, y, phi, levels=levels, linestyles='dashed',alpha=0.5)
        mean = np.mean(psi)
        std = np.std(psi)
        levels = np.linspace(norm(mean,std).ppf(0.1), norm(mean,std).ppf(0.9), num)  
        ax.contour(x, y, psi, levels=levels)
    else:
        if pot == True:
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
        
    if stream != False:
        zero = np.argwhere(np.around(psi,1) == stream)
        zeros = np.zeros((2,zero.shape[0]))
        for i in range(zero.shape[0]):
            zeros[0,i] = x[zero[i,0],zero[i,1]]
            zeros[1,i] = y[zero[i,0],zero[i,1]]
        ax.plot(zeros[0],zeros[1],'c*')
        
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
            
    fig1, ax1 = plt.subplots()
    ax1.contourf(x, y, mag)
            
if __name__ == "__main__":
    #HW4('Problem 1',func='lmao',num=50,pot=False,num2=40)
    #HW4('Problem 2',func='lmao',num=50,pot=False,num2=40)
    HW4('Problem 3',func='lmao',num=50,pot=False,num2=40)
    #HW4(func='lmao',num2=40)


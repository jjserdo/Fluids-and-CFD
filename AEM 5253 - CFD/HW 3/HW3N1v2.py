import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la

x0 = 5
D = 0.1

def Canal(x,t):
    return 1/(np.sqrt(4*np.pi*D*t)) * np.exp(-(x-x0)**2/(4*D*t))

N = 200

x = np.linspace(0,10,N+1)
dx = 10/N
#print(dx)
dt = np.round(dx**2/(2*D),10)
#dt = 0.0142 # hmm only stops working at 0.0142
print(dt)
numt = int(np.floor(2/dt))
#print(numt)
t = np.linspace(2,4,numt+1)
#print(t[0])
#print(t[1])
print('CFL is', str(4*dt/dx))

def bound(f):
    f[0] = 0
    f[-1] = 0
    return f

cti = Canal(x,2)
bound(cti)
ctf = Canal(x,4)
bound(ctf)

''' Right Hand Side Functions '''
def CS(f, n, i):
    RHS = D * (f[n][i+1] - 2*f[n][i] + f[n][i-1] ) / (dx)**2
    return RHS
''' Left Hand Side Functions '''
def Rk4(f, n, i):
    k0 = CS(f,n,i)
    k1 = CS(f,n,i) + dt/2 * k0
    k2 = CS(f,n,i) + dt/2 * k1 
    k3 = CS(f,n,i) + dt   * k2 
    df  = (k0 + 2*k1 + 2*k2 + k3) / 6
    LHS = f[n][i] + dt * df
    return LHS


# Exact Solution at t=2 and t=4
fig, ax = plt.subplots()
ax.plot(x, cti)
ax.plot(x, ctf)
plt.show()

C = np.zeros((len(t),len(x)))
C[0] = cti
for n in range(len(t)-1):
    for i in np.arange(N):
        C[n+1][i] = Rk4(C, n, i)
    C[n+1] = bound(C[n+1])
    
# Analytical Solution
plt.plot(x, C[0],'g',label='Analytical t=2')
#ax.plot(x, C[40],'g')
#ax.plot(x, C[80],'g')
#ax.plot(x, C[120],'g')
plt.plot(x, ctf, label='Analytical t=4')
plt.plot(x, C[-1],'r*', label='Rk4 at t=4')
plt.legend()
plt.title('Diffusion using RK4')
plt.xlabel('x')
plt.ylabel('Concetration, C')
plt.show()


#### Part 2 Crank Nicolson
dt = 1
# dt works even as high as 2 
numt = int(np.floor(2/dt))
print(numt)
t = np.linspace(2,4,numt+1)
print(t)
print('CFL is', str(4*dt/dx))

C = np.zeros((len(t),len(x)))
C[0] = cti
for n in np.arange(1,len(t)):#range(len(t)-1):
    A = np.zeros((len(x),len(x)))
    B = np.zeros(len(x))
    d = D*dt/(2*dx**2)
    B[0] = 0
    B[-1] = 0
    A[0][0] = 1
    A[-1][-1] = 1
    #print(C[0])
    for i in np.arange(1,len(B)-1):
        A[i,i-1] = -d
        A[i,i] = 1+2*d
        A[i,i+1] = -d
        B[i] = d*C[n-1][i-1] + (1-2*d)*C[n-1][i] + d*C[n-1][i+1]
    C[n] = la.solve(A,B)


      
# Analytical Solution
plt.plot(x, C[0],'g',label='Analytical t=2')
plt.plot(x, C[1], label = 'CN at t=3')
plt.plot(x, ctf, label='Analytical t=4')
plt.plot(x, C[-1],'r*', label='CN at t=4')
plt.legend()
plt.title('Diffusion using Crank-Nicolson using dt = 1')
plt.xlabel('x')
plt.ylabel('Concetration, C')
plt.show()







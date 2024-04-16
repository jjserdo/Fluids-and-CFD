"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1
"""
import numpy as np
import matplotlib.pyplot as plt
import math 

### Number 1 ###
ti = 0
tf = 15
yi = 4

# g(y(t), y)
def g(y, t):
    return - 2 * y;

# Analytical Solution solution
t_anal = np.linspace(ti,tf,1000)
print(t_anal)
y_anal = yi * np.exp(-2*t_anal)

# Time arrays
N01 = int((tf-ti)/0.1)
N05 = int((tf-ti)/0.5)
N10 = int((tf-ti)/1.0)
t_01 = np.linspace(ti,tf,N01+1)
t_05 = np.linspace(ti,tf,N05+1)
t_10 = np.linspace(ti,tf,N10+1)

# Explicit Euler
t = t_10
dt = 1
y_ex10 = np.zeros(len(t))
y_ex10[0] = yi

for i in range(len(t)-1):
    y_ex10[i+1] = y_ex10[i] + 1 * g(y_ex10[i],t[i])

fig1, ax1 = plt.subplots()
ax1.plot(t_anal,y_anal,'-bo',label='Analytical Solution') 
ax1.plot(t_ex10,y_ex10,'-ro',label='Explicit Euler')
ax1.set_title('Explicit Euler with Analytical Solution')
ax1.xlabel('Time')
ax1.ylabel('y')
ax1.set_xlim([0,15])
ax1.set_ylim([0, 4])
plt.show()

# Implicit Euler
t = t_10
dt = 1
y_im10 = np.zeros(len(t))
y_im10[0] = yi

k1 = 0
for i in range(len(t)-1):
    y_im10[i+1] = (y_im10[i] + 1 * (g(y_im10[i],0) + 2*y_im10[i]))/(1+2*dt)

fig2, ax2 = plt.subplots()
ax2.plot(t_anal,y_anal,'-bo',label='Analytical Solution') 
ax2.plot(t,y_im10,'-ro',label='Implicit Euler')
ax2.set_title('Implicit Euler with Analytical Solution')
ax2.xlabel('Time')
ax2.ylabel('y')
ax2.set_xlim([0,15])
ax2.set_ylim([0, 4])
plt.show()
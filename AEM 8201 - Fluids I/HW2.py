# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:12:26 2022

@author: jjser
"""

import numpy as np
import matplotlib.pyplot as plt

x0 = [0.5, 1, 1.5]
y0 = [0.5, 1, 1.5]

def stream(x0, y0):
    t = 3
    xf = x0 * (4+t) / 4
    x = np.linspace(x0, xf)
    y = y0 * x**2 / x0**2
    return x, y

fig, ax = plt.subplots()
for X0 in x0:
    for Y0 in y0:
        x, y = stream(X0, Y0)
        ax.plot(x,y,'-')
        
ax.set_title('Streamlines')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim([0,3])
ax.set_ylim([0,5])
fig.set_size_inches(5, 5)
plt.grid()
plt.savefig('images/hw2.png')
plt.show()
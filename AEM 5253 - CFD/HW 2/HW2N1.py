# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 17:03:28 2022

@author: jjser
"""
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,np.pi,100)
y = x
plt.plot(x,y)
y = np.sin(x)
plt.plot(x,y,label='S1')
y = np.sin(x)
plt.plot(x,y,label='S2')
y = (-np.sin(2*x)+8*np.sin(x))/6
plt.plot(x,y,label='S3')
y = (-np.sin(2*x)+8*np.sin(x))/6
plt.plot(x,y,label='S4')

title = "Real part of $k* \\Delta x$ vs $k \\Delta x$"
plt.title(title)
legend = plt.legend(loc='upper left',shadow=True,fontsize='large')
xlabel = "$k* \\Delta x$"
plt.xlabel(xlabel)
ylabel = "$k_{re} \\Delta x$"
plt.ylabel(ylabel)
plt.savefig('images/real.jpg')
plt.show()

y = np.zeros(len(x))
plt.plot(x,y)
y = -1+np.cos(x)
plt.plot(x,y,label='S1')
y = np.zeros(len(x))
plt.plot(x,y,label='S2')
y = (4*np.cos(x)-np.cos(2*x)-3)/6
plt.plot(x,y,label='S3')
y = np.zeros(len(x))
plt.plot(x,y,label='S4')

title = "Imaginary part of $k* \\Delta x$ vs $k \\Delta x$"
plt.title(title)
legend = plt.legend(loc='lower left',shadow=True,fontsize='large')
xlabel = "$k* \\Delta x$"
plt.xlabel(xlabel)
ylabel = "$k_{im} \\Delta x$"
plt.ylabel(ylabel)
plt.savefig('images/imag.jpg')
plt.show()

y = np.zeros(len(x)) + 10
plt.plot(x,y)

y = (x-np.sin(x)) / x * 100
plt.plot(x,y,label='S2')

y = (x-(-np.sin(2*x)+8*np.sin(x))/6) / x * 100
plt.plot(x,y,label='S4')

title = "Percent Error for changing $k \\Delta x$"
plt.title(title)
legend = plt.legend(loc='upper left',shadow=True,fontsize='large')
xlabel = "$k \\Delta x$"
plt.xlabel(xlabel)
ylabel = "% Error"
plt.ylabel(ylabel)
plt.savefig('images/error.jpg')
plt.show()


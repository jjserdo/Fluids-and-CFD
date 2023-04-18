import os
#path=os.getcwd()
#print(path)
os.chdir('C:/Users/jjser/OneDrive - Syracuse University/Fall 2022/AEM 5253/Project')

from MeshGenv7_plots import \
    plot_umag, plot_vel, get_actual \
    ,plot_density, plot_pressure, compute_cdrag
import pandas as pd
import numpy as np

u0 = 0.4
ga = 1.4
Uload = np.loadtxt('data/U_'+str(int(u0*100))+'.txt')
U = Uload.reshape(
    Uload.shape[0], Uload.shape[1] // 4, 4)
ruvp = U.copy()
ib = 1
jb = 1
ie = ib + U.shape[0] - 1 # ni
je = jb + U.shape[0] - 1 # nj
get_actual()
#print(U[0,:,1])
#print(ruvp[0,:,1])
plot_umag()


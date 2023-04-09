"""
Fall 2022 AEM 5253
Justine John "JJ" A. Serdoncillo
Homework 1 Number 5
"""
import numpy as np
from numpy import linalg as LA

Num = np.linspace(11,51,5)

def Agen(N):
    A = np.zeros((N,N))
    A[0, -1] = -1
    A[-1, 0] =  1
    for i in range(N):
        A[i-1,i  ] =  1
        A[i  ,i-1] = -1
    return A
def Bgen(N):
    B = np.zeros((N,N))
    B[0, -1] =  1
    B[-1, 0] =  1
    for i in range(N):
        B[i-1,i  ] =  1
        B[i  ,i-1] =  1
        B[i-1,i-1] =  2
    return B

for i in range(len(Num)):
    print('for N = ', str(Num[i]))
    Aeig = LA.eigvals(Agen(int(Num[i])))
    for j in range(4):
        AeigYeah = round(Aeig[j].real,2) + round(Aeig[j].imag,2) * 1j
        print(AeigYeah)
for i in range(len(Num)):
    print('for N = ', str(Num[i]))
    Beig = LA.eigvals(Bgen(int(Num[i])))
    Beig.argsort()[::-1]
    for j in range(4):
        BeigYeah = round(Beig[j].real,2) + round(Beig[j].imag,2) * 1j
        print(BeigYeah)



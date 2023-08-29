#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 15:48:22 2023

@author: homo
"""
import numpy as np
import cmath
from math import *
from matplotlib import pyplot as plt
from numpy import linalg as LA


#basic constants
a=1
b=1

p=2
q=100
fai=p/q

t_y = 1
t_x = 1



kx=0.1
ky=0.1

# matrix generator
def hamiltonian(kx,ky,q):
    ham=np.zeros((q,q),dtype=complex)
    fai = p/q
    for i in range(1,q-1):
        ham[i,i]=2*t_y*cos(2*pi*fai*i + ky*b)
        ham[i,i+1]=t_x
        ham[i,i-1]=t_x
        
    ham[0,1]=t_x
    ham[0,q-1]=cmath.exp(1j*q*a*kx)*t_x
    ham[0,0]=2*t_y*cos(ky*b)
    ham[q-1,0]=t_x*cmath.exp(-1j*q*a*kx)
    ham[q-1,q-2]=t_x
    ham[q-1,q-1]=2*t_y*cos(2*pi*fai*(q-1) + ky*b)
    
    return ham


#main part
def gcd(p,q):
    if q==0 :return p
    return gcd(q,p%q)

num =101
for p in range(1,num):
    for q in range(1,num):
        if q>p:
            if gcd(p,q)==1:
                phil=np.zeros(q)
                phil[:]=p/q
                eigs,vec = LA.eigh(hamiltonian(kx,ky,q))
                plt.plot(eigs,phil,'ko',markersize=0.2)

plt.ylabel('Flux (p/q)')    
plt.xlabel('Energy')
plt.title('Hofstadter spectrum')
plt.show()           
    



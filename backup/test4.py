# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, collimator, wedge_X, sec_bend, bending, drift2, quad2, sec_bend2, bending2,collimator2





nb_part=10
nb_pts_z = 100

beam = np.empty((nb_pts_z,7,nb_part)) 

refE = 160
Brho  = 1/300*sqrt(refE**2+2*938*refE)
   

for i in range(0,nb_part):
    
    z=0  
    
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    #E = 160 + np.random.normal(0,5)
    E = 160 + 3*np.random.choice([-1,0,1])
    
    #divX = 0.05
    #divY = 0.05
    divX=0
    divY=0
    E = 160
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,E])
    
    it_z = 0
    
    #ref_x[it_z] = 0
    
    
    
    L = 0.1
    N_segments = 1
    [beam[:,:,i],it_z] = drift2(L,beam[:,:,i],it_z,N_segments)
    
    
    
    
    
    [beam[:,:,i],it_z] = collimator2(0.1,.0,'lexan',beam[:,:,i],it_z,50)
    
    
    it_z += 1
    
    

plt.figure(0)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')

plt.figure(2)
plt.hist(beam[0,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_in")
plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')

plt.figure(3)
plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')
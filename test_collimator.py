# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, sec_bend, bending, collimator





nb_part=1000
nb_pts_z = 1000

beam = np.empty((nb_pts_z,7,nb_part)) 



refE = 160
Brho  = 1/300*sqrt(refE**2+2*938*refE)
   

for i in range(0,nb_part):
    
    z=0  
    
    
    divX = np.random.normal(0,0.02)
    divY = np.random.normal(0,0.02)
    #E = 160 + np.random.normal(0,5)
    E = 160 + 3*np.random.choice([-1,0,1])
    
    #divX = 0.05
    #divY = 0.05
    
    x = np.random.uniform(-0.002,0.002)
    y = np.random.uniform(-0.002,0.002)
    
    x = np.random.choice([-0.01,0,0.01])
    
    #x=0
    y=0
    divX=0
    divY=0
    E = 160
    beam[0,:,i] = np.array([z,x,divX,y,divY,0,E])
    
    it_z = 0
    
    #ref_x[it_z] = 0
    
    
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    
    N_segments = 100
    [beam[:,:,i],it_z] = collimator(0.16,0.001,'lexan',beam[:,:,i],it_z,N_segments)
    
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
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
#plt.hist(beam[0,6,:],20,alpha=0.5,label="E_in")
plt.hist(beam[it_z-1,6,:],20,alpha=0.5,label="E_out")
plt.legend(loc='upper right')


plt.figure(3)
plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')


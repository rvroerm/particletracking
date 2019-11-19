# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math as math
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, sec_bend, bending, wedge_X, solenoid, collimator





nb_part=1000
nb_pts_z = 1000

beam = np.empty((nb_pts_z,7,nb_part)) 

refE = 160
Brho  = 1/300*sqrt(refE**2+2*938*refE)
   

for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    divX = 0.05*np.random.choice([-1,0,1])
    
    #E = 160 + np.random.normal(0,3)
    E = 160 + 3*np.random.choice([-1,0,1])
    
    #divX = 0.05
    #divY = 0.05
    #E = 160
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,E])
    
    it_z = 0
    
    #ref_x[it_z] = 0
    
    
    # filter directly partcles above 50mrad
    L = 0.1
    [beam[:,:,i],it_z] = collimator(L,0.005*sqrt(2),'tantalum',beam[:,:,i],it_z,10)
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    L=0.15
    B=-0.96
    a=0.015
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    L=0.33
    B=1.03
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    L=0.21
    B=-0.9
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    L=0.1
    B=1
    pole_in = 0
    pole_out = 0
    gap = 0.01
    N_segments = 10
    #[beam[:,:,i],it_z] = sec_bend(L,B,0,beam[:,:,i],160,True,it_z,N_segments)
    [beam[:,:,i],it_z] = bending(L,B,0,pole_in,pole_out,gap,0.5,0,beam[:,:,i],160,True,it_z,N_segments)
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    use_wedge = True
    L = 0.01
    N_segments = 5
    if use_wedge: 
        [beam[:,:,i],it_z] = wedge_X(2*L,0.1,'lexan',beam[:,:,i],it_z,N_segments)
    else:
        [beam[:,:,i],it_z] = drift(2*L,beam[:,:,i],it_z,N_segments)
    
    L = 0.2
    [beam[:,:,i],it_z] = collimator(L,0.1*sqrt(2),'tantalum',beam[:,:,i],it_z,10)
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    new_refE = 154.17
    
    L=0.2
    B=-1.5
    pole_in = 60
    pole_out = 60
    gap = 0.01
    N_segments = 50
    #[beam[:,:,i],it_z] = sec_bend(L,B,0,beam[:,:,i],160,True,it_z,N_segments)
    [beam[:,:,i],it_z] = bending(L,B,0,pole_in,pole_out,gap,0.5,0,beam[:,:,i],new_refE,True,it_z,N_segments)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    L=0.2
    B=-1.
    a=0.1
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
#    
#    L=0.7
#    B=10
#    Brho  = 1/300*math.sqrt(E**2+2*938*E)
#    if Brho/B < L*10 :
#        # need more segments
#        N_segments = max(int(L/(Brho/B))*10,1)
#    else:
#        N_segments = 1
#    
#    [beam[:,:,i],it_z] = solenoid(L,B,0.1,beam[:,:,i],it_z,N_segments)
#    
#    
#    
#    L = 0.5
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    
    
    
    
    
    
    
    
    
    
#    L=0.4
#    B=1
#    a=0.1
#    N_segments = 10
#    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
#    
#    L = 0.3
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
#    
#    
#    L=0.3
#    B=-1
#    a=0.1
#    N_segments = 10
#    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
#    
#    
#    L = 0.2
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
#    
#    
#    L=0.2
#    B=1
#    a=0.1
#    N_segments = 10
#    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments)
#    
#    L = 1
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
#    
#    L = 1
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)


    it_z += 1




plt.figure(0)
plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')

plt.figure(10)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,2,:])
plt.title('divX')
plt.grid(which='major')



plt.figure(1)
plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')

plt.figure(11)
plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')


print("E mean = ",np.mean(beam[it_z-1,6,:]))

plt.figure(2)
plt.hist(beam[0,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_in")
plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')



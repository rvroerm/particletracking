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
from transfer_functions import drift, quad, sec_bend, bending, wedge_X, solenoid, collimator,EtoP,PtoE,sextupole


plt.close('all')


nb_part=100
nb_pts_z = 1000

beam = np.empty((nb_pts_z,7,nb_part)) 

refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)
   

for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    divX = 0.05*np.random.choice([-1,0,1])
    
    #divX=0.05
    divY=0.05
    
    
    
    E = 160 + 30*np.random.choice([-1,0,1])
    #E = 130 + 30*i
    E = 160 + np.random.normal(0,3)
    E = 160 + np.random.uniform(-30,30)
    #E = 160
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    
    divX = 0.000001
    divY = 0.000001
    
    beam[0,:,i] = np.array([z,.01,divX,.01,divY,0,dponp])
    
    it_z = 0
    
    #ref_x[it_z] = 0
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,10)
    
    L=0.2
    B=0.7
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,refE,N_segments)
    
    
    L = 0.4
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,10)
    
    
    L = 0.1
    B=0.8
    a=0.02
    alpha = 30
    N_segments = 10
    use_sext = True    
    if use_sext :
        [beam[:,:,i],it_z] = sextupole(L,B,a,alpha,beam[:,:,i],it_z,refE,N_segments)
    else:
        [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
#    
#    L=0.1
#    B = 1
#    a=0.04
#    N_segments = 10
#    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,refE,N_segments)
#    
    
    L = 0.5
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,10)
    
#plt.figure(23)
#plt.plot(beam[it_z,1,:]*1000,beam[it_z,2,:]*1000,marker='.', linestyle='None',label="emittance X")
#plt.plot(beam[it_z,3,:]*1000,beam[it_z,4,:]*1000,marker='.', linestyle='None',label="emittance Y")
#plt.legend(loc='upper right')
#
#
plt.figure(0)
#plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')

#plt.figure(10)
#plt.plot(beam[0:it_z,0,:],beam[0:it_z,2,:])
#plt.title('divX')
#plt.grid(which='major')



plt.figure(1)
#plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')


#plt.figure(25)
#plt.hist(beam[it_z-1,5,:],20,alpha=0.5,label="dL")
##plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
#plt.legend(loc='upper right')
#
##plt.figure(11)
##plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
##plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
###plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
##plt.legend(loc='upper right')
#
#
#
#
#plt.figure(2)
#plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-20,refE+10),alpha=0.5,label="E_in")
#plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,:])),100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
#plt.legend(loc='upper right')
#
#
#print(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,:])))
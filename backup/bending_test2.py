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
from transfer_functions import drift, quad, sec_bend, bending, wedge_X, solenoid, collimator,EtoP,PtoE


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
    
    
    
    E = 160 + 3*np.random.choice([-1,0,1])
    #E = 160 + np.random.normal(0,3)
    #E = 160 + np.random.uniform(-30,30)
    #E = 157
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    
    #divX = 0.05
    #divY = 0.05
    
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    it_z = 0
    
    #ref_x[it_z] = 0
    
    
    # filter directly partcles above 50mrad
    L = 0.1
    [beam[:,:,i],it_z] = collimator(L,0.005*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    L=0.15
    B=-0.96
    a=0.015
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=0.95
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.25
    B=-0.7
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    L=0.7
    B=1.5
    pole_in = 0
    pole_out = 0
    gap = 0.01
    N_segments = 10
    #[beam[:,:,i],it_z] = sec_bend(L,B,0,beam[:,:,i],160,True,it_z,N_segments)
    [beam[:,:,i],it_z] = bending(L,B,0,pole_in,pole_out,gap,0.5,0,beam[:,:,i],refE,it_z,N_segments)
    
    L = 0.24
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)




last_itz = it_z 
plt.figure(13)
plt.plot(beam[it_z,1,:]*1000,beam[it_z-1,2,:]*1000,marker='.', linestyle='None',label="emittance X")
#plt.plot(beam[it_z,3,:]*1000,beam[it_z-1,4,:]*1000,marker='.', linestyle='None',label="emittance Y")
plt.legend(loc='upper right')
   
    
for i in range(0,nb_part):
    
    it_z = last_itz 
        
    
    use_wedge = True
    L = 0.008
    N_segments = 5
    if use_wedge: 
        [beam[:,:,i],it_z] = wedge_X(2*L,0.0098,'lexan',beam[:,:,i],it_z,N_segments,refE)
    else:
        [beam[:,:,i],it_z] = drift(2*L,beam[:,:,i],refE,it_z,N_segments)
        
    
    
    
    
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
#cut to see distributions
new_refE = np.mean(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
print("E mean = ",new_refE)
print("p mean = ",np.mean((ref_p*(1+beam[it_z,6,:]))))


plt.figure(11)
plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')

plt.figure(12)
plt.hist(beam[it_z-1,2,:],20,alpha=0.5,label="divX")
plt.hist(beam[it_z-1,4,:],20,alpha=0.5,label="divY")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')

plt.figure(13)
plt.plot(beam[it_z,1,:]*1000,beam[it_z,2,:]*1000,marker='.', linestyle='None',label="emittance X out")
#plt.plot(beam[it_z,3,:]*1000,beam[it_z,4,:]*1000,marker='.', linestyle='None',label="emittance Y out")
plt.legend(loc='upper right')
plt.xlabel('size [mm]')  
plt.ylabel('divergence [mrad]')


last_itz = it_z    
    
for i in range(0,nb_part):
    
    it_z = last_itz 
    
    L = 0.1
    [beam[:,:,i],it_z] = collimator(L,0.025*sqrt(2),'tantalum',beam[:,:,i],it_z,10,new_refE)
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    L=0.15
    B=-0.5413
    a=0.015
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,new_refE)
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=0.976
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,new_refE)
    
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=-0.552
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,new_refE)
    
    
    L = 0.7
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=0.45861
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,new_refE)
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=-0.259
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,new_refE)
    
    L = 1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)


    #it_z += 1



plt.figure(23)
plt.plot(beam[it_z,1,:]*1000,beam[it_z,2,:]*1000,marker='.', linestyle='None',label="emittance X")
plt.plot(beam[it_z,3,:]*1000,beam[it_z,4,:]*1000,marker='.', linestyle='None',label="emittance Y")
plt.legend(loc='upper right')


plt.figure(0)
plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')


plt.figure(10)
plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,2,:])
plt.title('divX')
plt.grid(which='major')

#plt.figure(10)
#plt.plot(beam[0:it_z,0,:],beam[0:it_z,2,:])
#plt.title('divX')
#plt.grid(which='major')



plt.figure(1)
plt.ylim(-.1, 0.1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')


plt.figure(25)
plt.hist(beam[it_z-1,5,:],20,alpha=0.5,label="dL")
#plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')

#plt.figure(11)
#plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
#plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
##plt.hist(beam[it_z-1,6,:],100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
#plt.legend(loc='upper right')




plt.figure(2)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-20,refE+10),alpha=0.5,label="E_in")
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,:])),100,range=(refE-20,refE+10),alpha=0.5,label="E_out")
plt.legend(loc='upper right')
#
#
#print(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,:])))
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import solenoid, drift,EtoP, PtoE, collimator


plt.close('all')

#print(math.pi)
a=np.arange(0,100,.01)
b=np.random.normal(0,a)
#plt.plot(a,b)
#plt.close
#print(np.size(a))

#b=np.random.normal(0,0.1,np.size(a))
#plt.plot(a,b)
#plt.hist(b)


####

#x=np.arange(0,10,.01)
#y=np.arange(0,1,.01)
#z=np.arange(0,10,.01)
#
#mat=[x,y,z]
#
#plt.plot(mat[0][0:100],mat[1][0:100])
#print(mat[1][5])


nb_part=1000
nb_pts_z=1000 #100 is nb of points vs z, to refine
#beam = np.empty(shape=[nb_part,100,6]) #100 is nb of points vs z, to refine

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

refE=160
ref_p = EtoP(refE)


for i in range(0,nb_part):
    it_z = 0
    
    z=0
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    divX = 0.1 + np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.1 + 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    #divX = 0.05
    
    E = 150
    E = np.random.uniform(E-20,E+20)
    #E = 160 + 15*np.random.choice([-1,0,1])
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    
    beam[it_z,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    
    L = 0.3
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    
    #it_z = it_z + 1
    L=0.8
    B=7
    Brho  = 1/300*sqrt(E**2+2*938*E)
    
    if Brho/B < L*10 :
        # need more segments
        N_segments = max(int(L/(Brho/B))*10,1)
    else:
        N_segments = 1
    
    [beam[:,:,i],it_z] = solenoid(L,B,0.1,beam[:,:,i],it_z,refE,N_segments)
    
    
    L = 1.25
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = 0.01
    [beam[:,:,i],it_z] = collimator(L,0.002*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    
    L = 0.1
    [beam[:,:,i],it_z] = collimator(L,0.005*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    
    L = 0.25
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    
    L=0.6
    B=7
    Brho  = 1/300*sqrt(E**2+2*938*E)
    
    if Brho/B < L*10 :
        # need more segments
        N_segments = max(int(L/(Brho/B))*10,1)
    else:
        N_segments = 1
    
    [beam[:,:,i],it_z] = solenoid(L,B,0.1,beam[:,:,i],it_z,refE,N_segments)
    
    L = 1
    [beam[:,:,i],it_z] = collimator(L,0.05*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    print(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,i])))
    
    
    
    

plt.figure(0)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')



plt.figure(2)
plt.hist(beam[it_z-1,1,:],20,alpha=0.5,label="X")
plt.hist(beam[it_z-1,3,:],20,alpha=0.5,label="Y")
plt.legend(loc='upper right')

plt.figure(3)
plt.hist(beam[it_z-1,5,:],20,alpha=0.5,label="dL")
plt.legend(loc='upper right')

plt.figure(4)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-60,refE+20),alpha=0.5,label="E_in")
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z-1,6,:])),100,range=(refE-60,refE+20),alpha=0.5,label="E_out")
plt.legend(loc='upper right')
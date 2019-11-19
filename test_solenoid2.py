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
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.1 + 0.05*np.random.choice([-1,0,1])
    
    divX = 0.1 + 0.05*np.random.choice([-1,0,1])
    divY = 0.05*np.random.choice([-1,0,1])
    
    divX = 0.05
    divY = 0.05
    
    E = 150
    E = np.random.uniform(E-20,E+20)
    E = 70 + 10*np.random.choice([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
    E = 70 + 70*np.random.choice([0,1,2])
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    
    beam[it_z,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    
    L = 0.3
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    
    #it_z = it_z + 1
    L=0.9
    B=15
    Brho  = 1/300*sqrt(E**2+2*938*E)
    N_segments = 50
    [beam[:,:,i],it_z] = solenoid(L,B,0.1,beam[:,:,i],it_z,refE,N_segments)
    
    
    L = 0.5
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    L=0.3
    B=15
    N_segments = 50
    [beam[:,:,i],it_z] = solenoid(L,B,0.1,beam[:,:,i],it_z,refE,N_segments)
    
    L = .3
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = .3
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    

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
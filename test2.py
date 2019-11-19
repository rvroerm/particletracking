# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from transfer_functions import combined_magnet, drift, EtoP, PtoE, Brho_scaling, bending


plt.close('all')


nb_part=1



radius = 0.03
k1 = 0.5
k2 = 0

old_refE = 150 # default energy considered for fields bleow 
refE = 150
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)


Brho_factor = Brho_scaling(old_refE,refE)    




nb_pts_z = 1000




###############################################################################
# initial conditions

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

for i in range(0,nb_part):
    
    
    z=0  
    it_z = 0
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    sizeX = 0.
    sizeY = 0.
    
    #divX = np.random.normal(0,0.1/2.35) #100 mrad FWMH
    #divY = np.random.normal(0,0.1/2.35)
    
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    divX = 0.0 
    divY = 0.0
    
    #divX=0
    #divY=0
    
    E = refE + 10*np.random.normal(0,1)
    E = refE + 1.5*np.random.uniform(-1,1)
    #E = refE + 1.5*np.random.choice([-1,0,1])
    E = refE + 10*(i-nb_part/2)/nb_part*2
    #E = refE +  1.5*(i-1)
    E = refE - 5
    
    
#    sizeX = np.random.uniform(-0.004,0.004)
#    sizeY = np.random.uniform(-0.004,0.004)
#    divX = np.random.uniform(-0.01,0.01)
#    divY = np.random.uniform(-0.01,0.01)
#    #sizeX = 0.004
#    #sizeY = 0.004
#    #divX = 0.01
#    #divY = 0.01
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,sizeX,divX,sizeY,divY,0,dponp])
    
    
###############################################################################    
# extraction section
last_z = beam[it_z,0,:]    
last_itz = it_z 


for i in range(0,nb_part):  
    it_z = last_itz 
    #[beam,it_z] = transport_input('transport_file_ESS.txt',beam,refE,i,N_segments,gap,k1,k2,last_z,last_itz,Brho_factor,kill_lost_particles)
    
    L = 0.3
    B=1
    [beam[:,:,i],it_z] = bending(L,B,0,0,0,0.004,k1,k2,beam[:,:,i],refE,it_z,N_segments = 10,kill_lost_particles=False)









plt.figure(0)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.ylim((-0.05,0.05))
plt.xlabel('z [m]')
plt.ylabel('x [m]')
#plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.ylim((-0.05,0.05))
plt.xlabel('z [m]')
plt.ylabel('y [m]')
plt.grid(which='major')
#plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))

plt.figure(4)
plt.hist(beam[it_z,1,:],100,range=(-0.03,0.03),alpha=0.5,label="X")
plt.hist(beam[it_z,3,:],100,range=(-0.03,0.03),alpha=0.5,label="Y")
plt.legend(loc='upper right')

plt.figure(41)
plt.scatter(beam[it_z,1,:],beam[it_z,3,:])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.xlim((-0.05,0.05))
plt.ylim((-0.05,0.05))








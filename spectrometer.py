# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:25:04 2019

@author: HILAM
"""

from math import sqrt, pi
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from transfer_functions import bending, EtoP, RF_cavity, collimator, quad, drift, PtoE, PtoV, Gaussian_fit,Brho_scaling, slit


plt.close('all')


nb_part=10000



gap = 0.02 
k1 = 0.5
k2 = 0

old_refE = 160 # default energy considered for fields bleow 
refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)


Brho_factor = Brho_scaling(old_refE,refE)    



########################################
N_segments = 10
nb_pts_z = 1000

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 
    



for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    #sizeX = 0.00001
    #sizeY = 0.00001
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    #divX = 0.05 
    #divY = 0.05
    
    #divX=0
    #divY=0
    
    E = refE + 10*np.random.normal(0,1)
    E = refE + 3*np.random.uniform(-1,1)
    E = refE + 30*np.random.choice([-1,0,1])
    E = refE/2 + 50*(i-nb_part/2)/nb_part*2 + 30
    #E = refE +  3*(i-1)
    #E = refE
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,sizeX,divX,sizeY,divY,0,dponp])
    
    
    it_z = 0
    
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = 0.05
    N_segments = 10
    [beam[:,:,i],it_z] = slit('Y',L,-0.002,0.002,'tantalum',beam[:,:,i],it_z,refE,N_segments)
    
    L = 0.05
    N_segments = 10
    [beam[:,:,i],it_z] = slit('X',L,-0.001,0.001,'tantalum',beam[:,:,i],it_z,refE,N_segments)
    
    L=0.2
    B=1.5*Brho_factor
    n=0
    N_segments = 10
    pole_in = 0
    pole_out = 0
    [beam[:,:,i],it_z] = bending(L,B,n,pole_in,pole_out,gap,k1,k2,beam[:,:,i],refE,it_z,N_segments)
    
    # cleaning slit
    L = 0.05
    N_segments = 10
    [beam[:,:,i],it_z] = slit('Y',L,-0.006,0.006,'tantalum',beam[:,:,i],it_z,refE,N_segments)
    
    L = 0.35
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = 0.01
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
plt.figure(0)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
#plt.ylim((-0.05,0.05))
plt.xlabel('z [m]')
plt.ylabel('x [m]')
#plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
#plt.ylim((-0.05,0.05))
plt.xlabel('z [m]')
plt.ylabel('y [m]')
plt.grid(which='major')


print('efficiency = ',(nb_part-np.isnan(beam[it_z,6,:]).sum())/nb_part*100,' %')


plt.figure(2)
plt.scatter(beam[it_z,1,:],beam[it_z,3,:])
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.figure(21)
plt.scatter(beam[it_z,1,:],np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])))
plt.xlabel('x [m]')
plt.ylabel('E_in [ MeV]')

plt.figure(22)
plt.scatter(beam[it_z,1,:],np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
plt.xlabel('x [m]')
plt.ylabel('E_out [ MeV]')

plt.figure(23)
plt.scatter(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
plt.xlabel('E_in [ MeV]')
plt.ylabel('E_out [ MeV]')
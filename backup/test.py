# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP, drift, bending


plt.close('all')


nb_part=3

input_file = "C:/TRANS/for001.dat"

#input_file = "transport_input.txt"
#input_file = "HIL GTR B1G90 B1.5T largeEacceptance 3 quad types corrected Q2 variation study.txt"
#input_file = "test.txt"
#input_file = "achromat_90.txt"

N_segments = 10
nb_pts_z = transport_count_lines(input_file,N_segments)

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

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
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    divX = 0.045
    divY = 0.045
    
    E = refE + 10*np.random.normal(0,1)
    E = refE + 3*np.random.choice([-1,0,1])
    E = refE   + 3*(i-1)
    
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    it_z = 0
    
    beam = transport_input(input_file,beam,refE,i,N_segments)
    
#    
#    L = 0.2
#    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
#    
#    
#    L=0.1
#    B=1
#    pole_in = 0
#    pole_out = 0
#    gap = 0.01
#    N_segments = 10
#    #[beam[:,:,i],it_z] = sec_bend(L,B,0,beam[:,:,i],160,True,it_z,N_segments)
#    [beam[:,:,i],it_z] = bending(L,B,0,pole_in,pole_out,gap,0.5,0,beam[:,:,i],refE,it_z,N_segments)
#    
#    
#    L = 0.2
#    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)

plt.figure(0)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,3,:])
plt.title('Y')
plt.grid(which='major')


plt.figure(2)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,2,:])
plt.title('divX')
plt.grid(which='major')
plt.figure(3)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,4,:])
plt.title('divY')
plt.grid(which='major')

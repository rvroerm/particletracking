# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP


plt.close('all')


nb_part=1000
#input_file = "transport_input.txt"
input_file = "transport_input_ESS.txt"
input_file = "test.txt"
#input_file = "achromat_8quads.txt"

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
    
    divX = 0.05*np.random.choice([-1,0,1])
    divY = 0.05*np.random.choice([-1,0,1])
    E = 160 + 3*np.random.normal(0,1)
    E = 160 + 10*np.random.choice([-1,0,1])
    
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    #divX = 0.05
    #divY = 0.05
    
    #E = 160
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    it_z = 0
    
    
    beam = transport_input(input_file,beam,refE,i,N_segments)
    
    
    

plt.figure(0)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,3,:])
plt.title('Y')
plt.grid(which='major')

plt.figure(3)
plt.hist(beam[it_z-1,5,:],20,alpha=0.5,label="dL")
plt.legend(loc='upper right')


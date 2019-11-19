# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP, drift, bending, GTR_layout_from_transport, PtoE, collimator


plt.close('all')


nb_part=1000

input_file = "C:/TRANS/for001.dat"
input_file = "D:/beamline/transport/Transport code/inputs/HIL test 3T magnets double achromat.txt"



#input_file = "transport_input.txt"
#input_file = "HIL GTR B1G90 B1.5T largeEacceptance 3 quad types corrected Q2 variation study.txt"
#input_file = "test.txt"
#input_file = "achromat_90.txt"

gap = 0.03 
k1 = 0.7
k2 = 0


refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)

#######################################
# compute GTR drawing
nb_pts_z = transport_count_lines(input_file,1) 
test = np.zeros(shape=[nb_pts_z,2]) 
test = GTR_layout_from_transport(input_file,test,refE)
plt.figure(55)
plt.scatter(test[0:nb_pts_z-1,0],test[0:nb_pts_z-1,1])



########################################
N_segments = 10
nb_pts_z = transport_count_lines(input_file,N_segments) 

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 
    



for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
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
    #E = refE + 3*np.random.choice([-1,0,1])
    #E = refE +  3*(i-1)
    
    
    print(E)
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,.00001,divX,.00001,divY,0,dponp])
    
    [beam,it_z] = transport_input(input_file,beam,refE,i,N_segments,gap,k1,k2,0,0)
    
    
    
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
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.xlabel('z [m]')
plt.ylabel('x [m]')
plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.xlabel('z [m]')
plt.ylabel('y [m]')
plt.grid(which='major')
#plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))


plt.figure(2)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,2,:])
plt.title('divX')
plt.grid(which='major')
plt.figure(3)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,4,:])
plt.title('divY')
plt.grid(which='major')

plt.figure(4)
plt.hist(beam[it_z,1,:],100,range=(-0.03,0.03),alpha=0.5,label="X")
plt.hist(beam[it_z,3,:],100,range=(-0.03,0.03),alpha=0.5,label="Y")
plt.legend(loc='upper right')

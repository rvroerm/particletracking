# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP, drift, bending, GTR_layout_from_transport, PtoE, collimator, Brho_scaling


plt.close('all')


nb_part=1000

input_file = "C:/TRANS/for001.dat"
#input_file = "D:/beamline/transport/Transport code/inputs/HIL test 3T magnets double achromat.txt"
#input_file = "D:/beamline/transport/Transport code/inputs/HIL 1.5T magnets double achromat.txt"

#input_file = "D:/beamline/transport/Transport code/inputs/quads_Dansinger.txt"
#input_file = "D:/beamline/transport/Transport code/inputs/temp.txt"


#input_file = "transport_input.txt"
#input_file = "HIL GTR B1G90 B1.5T largeEacceptance 3 quad types corrected Q2 variation study.txt"
#input_file = "test.txt"
#input_file = "achromat_90.txt"

gap = 0.02 
k1 = 0.7
k2 = 0

old_refE = 230 # default energy considered for fields bleow 
refE = 230
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)


Brho_factor = Brho_scaling(old_refE,refE)    



#######################################
# compute GTR drawing
nb_pts_z = transport_count_lines(input_file,1) 
layout = np.zeros(shape=[nb_pts_z,2]) 
layout = GTR_layout_from_transport(input_file,layout,old_refE)
plt.figure(55)
plt.scatter(layout[0:nb_pts_z-1,0],layout[0:nb_pts_z-1,1])



########################################
N_segments = 10
nb_pts_z = transport_count_lines(input_file,N_segments) 

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 
    



for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    #sizeX = 0.00001
    #sizeY = 0.00001
    
    #divX = np.random.normal(0,0.1/2.35) #100 mrad FWMH
    #divY = np.random.normal(0,0.1/2.35)
    
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    #divX = 0.05 
    #divY = 0.05
    
    #divX=0
    #divY=0
    
    E = refE + 10*np.random.normal(0,1)
    E = refE + 1.5*np.random.uniform(-1,1)
    #E = refE + 1.5*np.random.choice([-1,0,1])
    #E = refE + 10*(i-nb_part/2)/nb_part*2
    #E = refE +  1.5*(i-1)
    #E = refE-1
    
    
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
    
    [beam,it_z] = transport_input(input_file,beam,refE,i,N_segments,gap,k1,k2,0,0,Brho_factor)
    
    
    
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
#    #[beam[:,:,i],it_z] = sec_bend(L,B,0,gap,beam[:,:,i],160,True,it_z,N_segments)
#    [beam[:,:,i],it_z] = bending(L,B,0,pole_in,pole_out,gap,0.5,0,beam[:,:,i],refE,it_z,N_segments)
#    
#    
#    L = 0.2
#    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)

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

plt.figure(5)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E source")
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_out GTR")
plt.legend(loc='upper right')
plt.xlabel('E [MeV]')
plt.ylabel('nb of protons (arb. units)')








E_list_in = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:]))
DeltaE = 3

nb_part_in_ESS = ((E_list_in < refE+DeltaE) & (E_list_in > refE-DeltaE)).sum()

E_list_out = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
E_list_out[np.isnan(E_list_out)] = 0 


nb_part_out_GTR = ((E_list_out < refE+DeltaE) & (E_list_out > refE-DeltaE)).sum()

print('total efficiency within E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_GTR/nb_part_in_ESS*100,'%')
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, collimator, wedge_X, sec_bend, bending, drift, quad, sec_bend, bending, EtoP, PtoE, Brho_scaling


plt.close('all')


nb_part=1000
#input_file = "transport_input.txt"
input_file = "transport_input_ESS.txt"
#input_file = "test.txt"


N_segments = 1000
nb_pts_z = 1000

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

old_refE = 160 # default energy considered for fields bleow 
refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)


Brho_factor = Brho_scaling(old_refE,refE)    

   

for i in range(0,nb_part):
    
    z=0  
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    #divX = 0.05
    #divY = 0.05
    
    
    E = refE  + 3*np.random.normal(0,2.5)
    E = refE + 30*(i-nb_part/2)/nb_part*2
    #E = refE + np.random.uniform(-5,5)
    #E = 160 + 3*np.random.choice([-1,0,1])
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    #divX = 0.0
    #divY = 0.0
    
    #E = refE  + 3
    
    
    p = EtoP(E)
    dponp = (p-ref_p)/ref_p
    beam[0,:,i] = np.array([z,sizeX,divX,sizeY,divY,0,dponp])
    
    
    
    it_z = 0
    
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    L_D=0.2
    L_Q=0.2
    B=2*Brho_factor
    a=0.05
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,-B,a,beam[:,:,i],it_z,refE,N_segments)
    [beam[:,:,i],it_z] = drift(L_D,beam[:,:,i],refE,it_z,1)
    [beam[:,:,i],it_z] = quad(L_Q,B,a,beam[:,:,i],it_z,refE,N_segments)
    
    
    L = 1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = 0.1
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
#plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])))

plt.figure(4)
plt.hist(beam[it_z,1,:],100,range=(-0.03,0.03),alpha=0.5,label="X")
plt.hist(beam[it_z,3,:],100,range=(-0.03,0.03),alpha=0.5,label="Y")
plt.legend(loc='upper right')

plt.figure(41)
plt.scatter(beam[it_z,1,:],beam[it_z,3,:])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
#plt.xlim((-0.05,0.05))
#plt.ylim((-0.05,0.05))

plt.figure(5)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E source")
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_out GTR")
plt.legend(loc='upper right')
plt.xlabel('E [MeV]')
plt.ylabel('nb of protons (arb. units)')








E_list_in = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:]))
DeltaE = 6

nb_part_in_ESS = ((E_list_in < refE+DeltaE) & (E_list_in > refE-DeltaE)).sum()

E_list_out = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
E_list_out[np.isnan(E_list_out)] = 0 


nb_part_out_GTR = ((E_list_out < refE+DeltaE) & (E_list_out > refE-DeltaE)).sum()

print('total efficiency within E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_GTR/nb_part_in_ESS*100,'%')
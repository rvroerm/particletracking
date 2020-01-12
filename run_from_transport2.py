# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP, PtoE, Brho_scaling, split_transport_file, gaussian

from plot_beam_results import plot_beam







nb_part=1000

input_file = "C:/TRANS/for001.dat"
#input_file = "D:/beamline/transport/Transport code/inputs/HIL test 3T magnets double achromat.txt"


split_transport_file(input_file)

paraxial_correction = False
dpzTolerance = 10**-4

gap = 0.03
k1 = 0.5
k2 = 0

old_refE = 160 # default energy considered for fields bleow 
refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)


Brho_factor = Brho_scaling(old_refE,refE)    

kill_lost_particles = True
gapX = gap # case of CCT magnets




########################################
N_segments = 10
nb_pts_z = transport_count_lines(input_file,N_segments) 

###############################################################################
# initial conditions

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

for i in range(0,nb_part):
    
    
    z=0  
    it_z = 0
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    #sizeX = 0.00001
    #sizeY = 0.00001
    
    
    
    divX = np.random.normal(0,0.05/2.35*2) #+-50 mrad FWMH
    divY = np.random.normal(0,0.05/2.35*2)
    
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    #divX = 0.05*np.random.choice([-1,0,1])
    #divY = 0.05*np.random.choice([-1,0,1])
    
    #divX = 0.05 
    #divY = 0.05
    
    #divX=0
    #divY=0
    
    E = refE + 5*np.random.normal(0,1)
    E = refE + 1.5*np.random.uniform(-1,1)
    #E = refE + 1.5*np.random.choice([-1,0,1])
    E = refE + 5*(i-nb_part/2)/nb_part*2
    #E = refE +  1.5*(i-1)
    #E = refE
    
    
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
    
    [beam,it_z] = transport_input('transport_file_ESS.txt',beam,refE,i,N_segments,gap,k1,k2,last_z,last_itz,Brho_factor,kill_lost_particles,gap_X=gapX,paraxial_correction=paraxial_correction,dpz_tolerance=dpzTolerance)
    


E_list_in = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:]))
DeltaE = refE/100
nb_part_in_ESS = ((E_list_in < refE+DeltaE) & (E_list_in > refE-DeltaE)).sum()
nb_part_in_ESS2 = ((E_list_in < refE+2*DeltaE) & (E_list_in > refE-2*DeltaE)).sum()

E_list_out_ESS = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
E_list_out_ESS[np.isnan(E_list_out_ESS)] = 0 




nb_part_out_ESS = ((E_list_out_ESS < refE+DeltaE) & (E_list_out_ESS > refE-DeltaE)).sum()

print('ESS efficiency within E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_ESS/max(nb_part_in_ESS,1)*100,'%')


###############################################################################    
# gantry section
last_z = beam[it_z,0,:]
last_itz = it_z 
it_z_GTR = it_z

for i in range(0,nb_part):  
    it_z = last_itz 
    [beam,it_z] = transport_input('transport_file_GTR.txt',beam,refE,i,N_segments,gap,k1,k2,last_z,last_itz,Brho_factor,kill_lost_particles,gap_X=gapX,dpz_tolerance=dpzTolerance)

        








# make plots

plot_beam(input_file,beam,it_z,it_z_GTR,ref_p)




# output results

E_list_out_GTR = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
E_list_out_GTR[np.isnan(E_list_out_GTR)] = 0 


nb_part_out_GTR = ((E_list_out_GTR < refE+DeltaE) & (E_list_out_GTR > refE-DeltaE)).sum()

print('GTR efficiency within 1% E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_GTR/max(nb_part_out_ESS,1)*100,'%')
print('total efficiency within 1% E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_GTR/max(nb_part_in_ESS,1)*100,'%')

nb_part_out_GTR2 = ((E_list_out_GTR < refE+2*DeltaE) & (E_list_out_GTR > refE-2*DeltaE)).sum()
print('total efficiency within 2% E range [',refE-2*DeltaE,',',refE+2*DeltaE,']  = ',nb_part_out_GTR2/max(nb_part_in_ESS2,1)*100,'%')
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt, pi
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from transfer_functions import transport_input, transport_count_lines, EtoP, RF_cavity, collimator, quad, drift, PtoE, PtoV, Gaussian_fit



plt.close('all')


nb_part=1000
#input_file = "transport_input.txt"
input_file = "transport_input_ESS.txt"
#input_file = "test.txt"


N_segments = 1000
nb_pts_z = 1000

beam = np.empty(shape=[nb_pts_z,7,nb_part]) 

refE = 160
ref_p = EtoP(refE)
Brho  = 1/300*sqrt(refE**2+2*938*refE)
    

for i in range(0,nb_part):
    
    # initial conditions
    z=0  
    
    sizeX = np.random.uniform(-0.00001,0.00001)
    sizeY = np.random.uniform(-0.00001,0.00001)
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    divX = np.random.uniform(-0.05,0.05)
    divY = np.random.uniform(-0.05,0.05)
    
    E = refE  + 3*np.random.normal(0,2.5)
    E = refE + 10*(i-nb_part/2)/nb_part*2
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
    
    # filter directly partcles above 50mrad
    L = 0.11
    [beam[:,:,i],it_z] = collimator(L,0.05*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],it_z,1)
    
    
    
    L=0.15
    B=-0.967
    a=0.015
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    L=0.33
    B=0.99
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    
    L = 0.2
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.33
    B=-0.583
    a=0.04
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    
    L = 0.7
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.25
    B=0.2978
    a=0.02
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    L = 0.15
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L=0.25
    B=0.13472
    a=0.02
    N_segments = 10
    [beam[:,:,i],it_z] = quad(L,B,a,beam[:,:,i],it_z,N_segments,refE)
    
    L = 1.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    
    #L = 0.2
    #[beam[:,:,i],it_z] = collimator(L,0.005*sqrt(2),'tantalum',beam[:,:,i],it_z,10,refE)
    #[beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)




# RF buncher

last_itz = it_z 

z = beam[it_z,0,0]
dist_cavities = 0.07
cav_length = 0.005
aperture = 0.01/2
delta_E = 0.25

speed = PtoV(ref_p,938)

nb_cavities = 5
freq = speed/(dist_cavities+cav_length) 
phi = -2*pi*freq * z / speed + pi/2



print('freq = ',freq/10**6,' MHz') 
    
for i in range(0,nb_part):
    
    it_z = last_itz 
    
    
    
    #print(beam[it_z,:,i])
    #[beam[:,:,i],it_z] = drift(dist_cavities*nb_cavities,beam[:,:,i],refE,it_z,1)
    #print(beam[it_z,:,i])
    
    for it_cav in range(0,nb_cavities):
        #print(beam[it_z,:,i])
        [beam[:,:,i],it_z] = RF_cavity(delta_E,cav_length,freq,phi,beam[:,:,i],refE,it_z,10)
        #print(beam[it_z,:,i])
        #[beam[:,:,i],it_z] = RF_cavity_old(delta_E,freq,phi,beam[:,:,i],refE,it_z)
        #[beam[:,:,i],it_z] = drift(dist_cavities,beam[:,:,i],refE,it_z,1)
        [beam[:,:,i],it_z] = collimator(dist_cavities,aperture,'tantalum',beam[:,:,i],it_z,10,refE)
        #print(beam[it_z,6,i])
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    
    





plt.figure(0)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.grid(which='major')

#plt.figure(3)
#plt.hist(beam[it_z,5,:],20,alpha=0.5,label="dL")
#plt.legend(loc='upper right')


plt.figure(4)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E source")
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_in GTR")
plt.legend(loc='upper right')

#plt.figure(41)
#plt.plot(beam[0:it_z,0,:],np.vectorize(PtoE)(ref_p*(1+beam[0:it_z,6,:])))
#plt.title('E')
#plt.grid(which='major')



# calculate efficiency for a given energy window
E_list_in = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:]))
E_list_out = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
E_list_out[np.isnan(E_list_out)] = 0 
DeltaE = 3

nb_part_in_ESS = ((E_list_in < refE+DeltaE) & (E_list_in > refE-DeltaE)).sum()
nb_part_out_ESS = ((E_list_out < refE+DeltaE) & (E_list_out > refE-DeltaE)).sum()


print(nb_part_in_ESS)
print(nb_part_out_ESS)

print('ESS efficiency within E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_ESS/nb_part_in_ESS*100,'%')


# consider particles within 1cm only enter the GTR (collimator)
mask = ~np.isnan(beam[it_z,1,:])
nb_part_in_GTR = ((beam[it_z,1,:][mask]< 0.01) & (beam[it_z,1,:][mask] > -0.01) & (beam[it_z,3,:][mask]< 0.01) & (beam[it_z,3,:][mask] > -0.01)).sum()
print("")



plt.figure(51)
plt.scatter(beam[0,1,:],beam[0,2,:],alpha=0.5,label="emittance X in")
plt.scatter(beam[it_z,1,:],beam[it_z,2,:],alpha=0.5,label="emittance X out")
plt.legend(loc='upper right')


plt.figure(52)
plt.scatter(beam[0,3,:],beam[0,4,:],alpha=0.5,label="emittance X in")
plt.scatter(beam[it_z,3,:],beam[it_z,4,:],alpha=0.5,label="emittance X out")
plt.legend(loc='upper right')

##################################################################




# fit

#print('X')
#Gaussian_fit(beam[it_z,1,:],[100, 0., 0.001])
#print('divX')
#Gaussian_fit(beam[it_z,2,:],[30., 0., 0.001])
#print('Y')
#Gaussian_fit(beam[it_z,3,:],[100., 0., 0.001])
#print('divY')
#Gaussian_fit(beam[it_z,4,:],[30., 0., 0.001])


##################################################################
# GTR

input_file = "D:/beamline/transport/Transport code/inputs/GTR 3T magnets double achromat after buncher.txt"

gap = 0.03 
k1 = 0.7
k2 = 0

last_itz = it_z

for i in range(0,nb_part):
    it_z = last_itz 
    [beam,it_z] = transport_input(input_file,beam,refE,i,N_segments,gap,k1,k2,beam[it_z,0,0],it_z)
    L = 0.1
    [beam[:,:,i],it_z] = collimator(L,0.02,'tantalum',beam[:,:,i],it_z,10,refE)
    L = 0.1
    [beam[:,:,i],it_z] = drift(L,beam[:,:,i],refE,it_z,1)
    


plt.figure(10)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,1,:])
plt.title('X')
plt.ylim((-0.05,0.05))
plt.grid(which='major')
plt.figure(11)
plt.plot(beam[0:it_z,0,:],beam[0:it_z,3,:])
plt.title('Y')
plt.ylim((-0.05,0.05))
plt.grid(which='major')

plt.figure(55)
plt.hist(beam[it_z,1,:],100,range=(-0.03,0.03),alpha=0.5,label="X")
plt.hist(beam[it_z,3,:],100,range=(-0.03,0.03),alpha=0.5,label="Y")
plt.legend(loc='upper right')


plt.figure(4)
plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_out GTR")
plt.legend(loc='upper right')
plt.xlabel('E [MeV]')
plt.ylabel('nb of protons (arb. units)')

#nb_part_in = ((E_list_in < refE+DeltaE) & (E_list_in > refE-DeltaE)).sum()

#nb_part_out_GTR = ((beam[it_z,1,:]**2 + beam[it_z,3,:]**2) < 0.01**2).sum()

#E_list_out = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
#E_list_out[np.isnan(E_list_out)] = 0 


mask = ~np.isnan(beam[it_z,1,:])
nb_part_out_GTR = ((beam[it_z,1,:][mask]< 0.01) & (beam[it_z,1,:][mask] > -0.01) & (beam[it_z,3,:][mask]< 0.01) & (beam[it_z,3,:][mask] > -0.01)).sum()


print("")


print('GTR efficiency = ',nb_part_out_GTR/nb_part_in_GTR*100,'%')
#print('total efficiency within E range = ',nb_part_out_GTR/nb_part_in_ESS*100,'%')
print('total efficiency within E range [',refE-DeltaE,',',refE+DeltaE,']  = ',nb_part_out_GTR/nb_part_in_ESS*100,'%')
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math as math
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, collimator, wedge_circular



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


nb_part=10
nb_pts_z=100 #100 is nb of points vs z, to refine
#beam = np.empty(shape=[nb_part,100,6]) #100 is nb of points vs z, to refine

beam = np.zeros(shape=[nb_pts_z,7,nb_part]) 



for i in range(0,nb_part):
    it_z = 0
    
    z=0
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    
    beam[it_z,:,i] = np.array([z,.00001,divX,.00001,divY,0,160])
    
    
    it_z = it_z + 1
    length=0.1
    z=z+length
    beam[it_z,:,i] = z
    beam[it_z,1:7,i] = drift(length,beam[it_z-1,1:7,i])
    
    
    use_wedge = True
    L = 0.01
    if use_wedge:
        it_z = it_z + 1
        L_min=0.0011
        L_max=L
        z=z+L_max
        beam[it_z,0,i] = z
        beam[it_z,1:7,i] = wedge_circular(L_min,.004,L_max,0.035,'tantalum',10,beam[it_z-1,1:7,i])
    else:
        it_z = it_z + 1
        length=L
        z=z+length
        beam[it_z,:,i] = z
        beam[it_z,1:7,i] = drift(length,beam[it_z-1,1:7,i])
    
    
    it_z = it_z + 1
    length=0.5
    z=z+length
    beam[it_z,0,i] = z
    beam[it_z,1:7,i] = drift(length,beam[it_z-1,1:7,i])
    


plt.figure(0)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,1,:])
plt.title('X')
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:nb_pts_z,0,:],beam[0:nb_pts_z,3,:])
plt.title('Y')
plt.grid(which='major')

print(beam[it_z,1,:].std())
print(beam[it_z,3,:].std())

print("P66 = ",np.percentile(beam[it_z,1,:],66))
print("P80 = ",np.percentile(beam[it_z,1,:],80))
print("P90 = ",np.percentile(beam[it_z,1,:],90))
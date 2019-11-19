# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math as math
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import drift, quad, collimator, wedge



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


iterations=10

beam = np.empty(shape=[iterations,7])



for i in range(0,iterations):
    z=0
    
    divX = np.random.normal(0,0.05)
    divY = np.random.normal(0,0.05)
    #beam = np.array([.00001,divX,.00001,divY,0,160])
    beam[i,:] = np.array([z,.00001,divX,.00001,divY,0,160])
    
    
    #beam = np.array([.00001,.05,.00001,.05,0,160])
    #beam = np.array([.01,.0,.01,.0,0,160])
    #x_vs_z = np.array([z,beam[0]][i])
    
    length=0.1
    z=z+length
    beam[i,0] = z
    beam[i,1:7] = drift(length,beam[i,1:7])
    #x_vs_z = np.vstack((x_vs_z,np.array([z,beam[0][i]])))
    
    #length=0.15
    #z=z+length
    #beam[0:6] = quad(length,-0.964,0.015,beam[0:6])
    #x_vs_z = np.vstack((x_vs_z,np.array([z,beam[0]])))
    
    L_min=0.001
    L_max=0.01
    z=z+L_max
    beam[i,0] = z
    beam[i,1:7] = wedge(L_min,.005,L_max,.02,'tantalum',10,beam[i,1:7])
    #x_vs_z = np.vstack((x_vs_z,np.array([z,beam[0][i]])))
    
    
# =============================================================================
#     L=0.0001
#     z=z+L
#     beam = collimator(L,.005,'tantalum',beam)
#     x_vs_z = np.vstack((x_vs_z,np.array([z,beam[0]])))
# =============================================================================
    
    
    length=0.1
    z=z+length
    beam[i,0] = z
    beam[i,1:7] = drift(length,beam[i,1:7])
    #x_vs_z = np.vstack((x_vs_z,np.array(beam[i,1:7])))
    
    #length=0.1
    #z=z+length
    #beam = drift(length,beam)
    #x_vs_z = np.vstack((x_vs_z,np.array([z,beam[0]])))
    
    plt.plot(beam[:,0],beam[:,1])
    #print(beam[5])
    #del x_vs_z


 
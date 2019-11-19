# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt, atan, cos, sin
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from transfer_functions import EOM_particle_in_Bfield, PtoV, EtoP

import pandas as pd
import time
start_time = time.clock()




plt.close('all')


path = "D:/beamline/magnets/field maps/test/"
field_map = pd.read_csv(path + "test.csv")
spatial_res = 10


# replacing blank spaces with '_'  
field_map.columns =[column.replace(" ", "_") for column in field_map.columns] 
field_map.columns =[column.replace("[", "") for column in field_map.columns] 
field_map.columns =[column.replace("]", "") for column in field_map.columns] 

#print(field_map[['z_[mm]']])

x = 0
y =0
z = 50


temp_df = field_map.loc[(field_map['x_mm'] >= x - spatial_res) & (field_map['x_mm'] <= x + spatial_res) & (field_map['y_mm'] >= y - spatial_res) & (field_map['y_mm'] <= y + spatial_res) & (field_map['z_mm'] >= z - spatial_res) & (field_map['z_mm'] <= z + spatial_res)]
temp_df['distXY'] = (temp_df['x_mm'] - x)**2 + (temp_df['y_mm'] - y)**2 
temp_df['dist'] = (temp_df['x_mm'] - x)**2 + (temp_df['y_mm'] - y)**2 + (temp_df['z_mm'] - z)**2 
temp_df = temp_df.sort_values(['dist','distXY']) 


[Bx,By,Bz] = temp_df.iloc[0].loc[['Bx_T','By_T','Bz_T']]





nb_part=1

nb_it_z = 300
Delta_z = 0.001


beam = np.empty(shape=[nb_it_z,nb_part,6]) 

beam_ref = np.empty(shape=[nb_it_z,6]) 
costheta_mat = np.empty(shape=[nb_it_z]) 
sintheta_mat = np.empty(shape=[nb_it_z]) 


refE = 150
refp = EtoP(refE)
ref_speed = PtoV(refp,938)


# trajectory of reference particle
for it_z in range(0,nb_it_z):
        
    Delta_t = Delta_z/ref_speed
    
    By = 1
    
    [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz] = EOM_particle_in_Bfield(Delta_t,0,0,0,0,0,ref_speed,0,By)
    beam_ref[it_z,:] = [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz]
    #theta_ZX = -atan(ref_vx/ref_vz)
    costheta_mat[it_z] = 1/sqrt(1+(ref_vx/ref_vz)**2)
    sintheta_mat[it_z] = ref_vx/ref_vz/sqrt(1+(ref_vx/ref_vz)**2)
    



for it_p in range(0,nb_part):
    
    x=0
    y=0
    z=0
    vx = 0
    vy = 0
    vz = 1.5*10**8
    v0 = vz
    beam[it_p,0,:] = [x,y,z,vx,vy,vz]
    
    
    
    ref_x = 0
    ref_y = 0
    ref_z = 0
    
    
    
    for it_z in range(0,nb_it_z):
        
        #print(it_z)
        #By = 1
        
        
        # extract subset around region of interest
        temp_df = field_map.loc[(field_map['x_mm'] >= x - spatial_res) & (field_map['x_mm'] <= x + spatial_res) & (field_map['y_mm'] >= y - spatial_res) & (field_map['y_mm'] <= y + spatial_res) & (field_map['z_mm'] >= z - spatial_res) & (field_map['z_mm'] <= z + spatial_res)]
        # compute distance to point of interest
        temp_df['distXY'] = (temp_df['x_mm'] - x)**2 + (temp_df['y_mm'] - y)**2 
        temp_df['dist'] = (temp_df['x_mm'] - x)**2 + (temp_df['y_mm'] - y)**2 + (temp_df['z_mm'] - z)**2 
        temp_df = temp_df.sort_values(['dist','distXY']) 
        # get field at closest point
        [Bx,By,Bz] = temp_df.iloc[0].loc[['Bx_T','By_T','Bz_T']]
        
        
        
        Delta_t = Delta_z/vz
        [x,y,z,vx,vy,vz] = EOM_particle_in_Bfield(Delta_t,x,y,z,vx,vy,vz,Bx,By,Bz)
        
        
        # compute difference with reference particle
        x = x - beam_ref[it_z,0]
        y = y - beam_ref[it_z,1]
        
        vx = costheta_mat[it_z]*vx - sintheta_mat[it_z]*vz
        vz = sintheta_mat[it_z]*vx + costheta_mat[it_z]*vz
        
        # normalize to correct rounding errors
        v = sqrt(vx**2 + vy**2 + vz**2)
        vx = vx * v0/v
        vy = vy * v0/v
        vz = vz * v0/v
        
        
        beam[it_z,it_p,:] = [x,y,z,vx,vy,vz]
    

print("--- execution time: %s seconds ---" % (time.clock() - start_time))



plt.figure(0)
plt.plot(beam[:,0,2],beam[:,:,0])
plt.figure(1)
plt.plot(beam[:,0,2],beam[:,:,1])

plt.figure(2)
plt.plot(beam[:,0,2],beam[:,:,3],label="vx")
plt.plot(beam[:,0,2],beam[:,:,4],label="vy")
plt.plot(beam[:,0,2],beam[:,:,5],label="vz")
plt.legend(loc='upper right')


plt.figure(3)
plt.plot(beam[:,:,3],label="vx")
plt.plot(beam[:,:,4],label="vy")
plt.plot(beam[:,:,5],label="vz")
plt.legend(loc='upper right')



plt.figure(5)


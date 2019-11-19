# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from math import sqrt
import math as math
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
from transfer_functions import EOM_particle_in_Bfield, PtoV, EtoP

import pandas as pd
import time
start_time = time.perf_counter()





def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]

def find_nearest_sorted(array,value):
    # finds nearest value in sorted array
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return [array[idx-1],idx-1]
    else:
        return [array[idx],idx]


plt.close('all')


path = "D:/beamline/magnets/field maps/test/"
field_map = pd.read_csv(path + "test.csv")
spatial_res = 10


# replacing blank spaces with '_'  
field_map.columns = [column.replace(" ", "_") for column in field_map.columns] 
field_map.columns = [column.replace("[", "") for column in field_map.columns] 
field_map.columns = [column.replace("]", "") for column in field_map.columns] 

#print(field_map[['z_[mm]']])

x = 0
y =0

field_map[['x_mm','y_mm','z_mm']] /= 1000 # convert mm to m
field_map.columns = [column.replace("_mm", "_m") for column in field_map.columns]


nb_part=10

nb_it_z = 300
Delta_z = 0.001


beam = np.empty(shape=[nb_it_z,nb_part,6]) 

beam_ref = np.empty(shape=[nb_it_z,6]) 
costheta_mat = np.empty(shape=[nb_it_z]) 
sintheta_mat = np.empty(shape=[nb_it_z]) 


refE = 150
refp = EtoP(refE)
ref_speed = PtoV(refp,938)


#List unique values in z
z_list = field_map['z_m'].unique()



# trajectory of reference particle
for it_z in range(0,nb_it_z):
        
    Delta_t = Delta_z/ref_speed
    
    By = 1
    
    [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz] = EOM_particle_in_Bfield(Delta_t,0,0,0,0,0,ref_speed,0,By)
    beam_ref[it_z,:] = [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz]
    #theta_ZX = -atan(ref_vx/ref_vz)
    costheta_mat[it_z] = 1/sqrt(1+(ref_vx/ref_vz)**2)
    sintheta_mat[it_z] = ref_vx/ref_vz/sqrt(1+(ref_vx/ref_vz)**2)
    






for it_z in range(0,nb_it_z):    
    
    z = Delta_z*it_z
    
    
    # 1) extract 2D map from 3D map by selecting good z value and convert to numpy matrix
    [z_map,ind_z] = find_nearest_sorted(z_list,z)
    field_mapZ = field_map.loc[field_map['z_m'] == z_map].to_numpy()
    
    
    # 2) find trajectory of reference particle
    
    # find closest point in XY map
    point_index = distance.cdist([np.array([0,0])], field_mapZ[:,0:2]).argmin()
    # extract field at that point
    Bx = field_mapZ[point_index,3]
    By = field_mapZ[point_index,4]
    Bz = field_mapZ[point_index,5]
    [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz] = EOM_particle_in_Bfield(Delta_t,0,0,0,0,0,ref_speed,0,By)
    beam_ref[it_z,:] = [ref_x,ref_y,ref_z,ref_vx,ref_vy,ref_vz]
    
    # compute rotation coefficients
    costheta_mat[it_z] = 1/sqrt(1+(ref_vx/ref_vz)**2)
    sintheta_mat[it_z] = ref_vx/ref_vz/sqrt(1+(ref_vx/ref_vz)**2)
    
    
    
    # 3) loop for other particles
    
    for it_p in range(0,nb_part):
        
        if it_z == 1:
            x=0
            y=0
            z=0
            vx = 0
            vy = 0
            vz = 1.5*10**8
            v0 = vz
            beam[it_p,0,:] = [x,y,z,vx,vy,vz]
            
        
        
        
        # find closest point in XY map
        point_index = distance.cdist([np.array([x,y])], field_mapZ[:,0:2]).argmin()
        # extract field at that point
        Bx = field_mapZ[point_index,3]
        By = field_mapZ[point_index,4]
        Bz = field_mapZ[point_index,5]
        
        
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
    

print("--- execution time: %s seconds ---" % (time.perf_counter() - start_time))



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

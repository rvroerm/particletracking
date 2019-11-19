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
from transfer_functions import EOM_particle_in_Bfield, PtoV, EtoP, drift

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

def open_B_map(path,filename):
    field_map = pd.read_csv(path + filename)
    
    # replacing blank spaces with '_'  
    field_map.columns = [column.replace(" ", "_") for column in field_map.columns] 
    field_map.columns = [column.replace("[", "") for column in field_map.columns] 
    field_map.columns = [column.replace("]", "") for column in field_map.columns] 
    field_map.columns = [column.replace("/", "_") for column in field_map.columns] 
    
    field_map[['x_mm','y_mm','z_mm']] /= 1000 # convert mm to m
    field_map.columns = [column.replace("_mm", "_m") for column in field_map.columns]
    
    return field_map

def particle_in_B_map(field_map,beam,Delta_z,itz0,nb_it_z,it_p,refp,resol_x,resol_y,spatial_res = 10):
    # to do: get resolution from map and get rid of resol_x, resol_y
    
    
    beam_ref = np.empty(shape=[nb_it_z,7]) 
    
    
    
    #List unique values in z
    z_list = field_map['z_m'].unique()
    
    z = beam[itz0,0]
    
#    costheta_mat = np.empty(shape=[nb_it_z]) 
#    sintheta_mat = np.empty(shape=[nb_it_z]) 
#    
    for it_z in range(itz0+1,itz0+nb_it_z):    
        
        z = z + Delta_z
        
        # 1) extract 2D map from 3D map by selecting good z value and convert to numpy matrix
        [z_map,ind_z] = find_nearest_sorted(z_list,z)
        field_mapZ = field_map.loc[field_map['z_m'] == z_map].to_numpy()
        
        
        # 2) find trajectory of reference particle 
        
        # find point near the center in XY map
        point_index = distance.cdist([np.array([0,0])], field_mapZ[:,0:2]).argmin()
        # extract field at that point
        Bx = field_mapZ[point_index,3] 
        By = field_mapZ[point_index,4] 
        Bz = field_mapZ[point_index,5]
        
        
            
        
        beam_ref[it_z-itz0,:] = EOM_particle_in_Bfield(Delta_z,[z,0,0,0,0,0,0],refp,Bx,By,Bz)
        
        
#        # compute rotation coefficients on X axis (case of dipole component)
#        costheta_mat[it_z-itz0] = 1/sqrt(1+beam_ref[it_z-itz0,2]**2)
#        sintheta_mat[it_z-itz0] = beam_ref[it_z-itz0,2]/sqrt(1+beam_ref[it_z-itz0,2]**2)
#        
        
        # 3) Actual particle
        
            
        if it_z == 0:
            print("beam should be initialized before fumction, itz0>0")
        else:
            x = beam[it_z-1,1]
            y = beam[it_z-1,3]
            
        
        # find closest point in XY map
        point_index = distance.cdist([np.array([x,y])], field_mapZ[:,0:2]).argmin()
        # extract field at that point
        x_map = field_mapZ[point_index,0]
        y_map = field_mapZ[point_index,1]
        Bx = field_mapZ[point_index,3] + field_mapZ[point_index,6] * (y-y_map)
        By = field_mapZ[point_index,4] + field_mapZ[point_index,7] * (x-x_map)
        Bz = field_mapZ[point_index,5]
        
        if np.isnan([Bx,By,Bz]).any():
            print('error: field is not valid')
            print('z = ',z)
            print('[x,y] = ',[x,y])
            print('[Bx,By,Bz] = ',[Bx,By,Bz])
            
        
        
        #v0 = sqrt(vx**2 + vy**2 + vz**2)
        
        beam[it_z,:] = EOM_particle_in_Bfield(Delta_z,beam[it_z-1,:],refp,Bx,By,Bz)
        
        # compute difference with reference particle
        beam[it_z,1] = beam[it_z,1] - beam_ref[it_z-itz0,1]
        beam[it_z,2] = beam[it_z,1] - beam_ref[it_z-itz0,2]
        
        #vx = costheta_mat[it_z-itz0]*vx - sintheta_mat[it_z-itz0]*vz
        #vz = sintheta_mat[it_z-itz0]*vx + costheta_mat[it_z-itz0]*vz
        
        ## normalize to correct rounding errors
        #v = sqrt(vx**2 + vy**2 + vz**2)
        #vx = vx * v0/v
        #vy = vy * v0/v
        #vz = vz * v0/v
        
        #beam[it_z,it_p,:] = [x,y,z,vx,vy,vz]
            
    return [beam,it_z]



plt.close('all')




nb_part=100

nb_it_z = 30
Delta_z = 0.01

nb_pts_z = 100

beam = np.empty(shape=[nb_pts_z,nb_part,7]) 

refE = 160

path = "D:/beamline/magnets/field maps/test/"
filename = "test_quad_B1T_aperture15_noise1e-4.csv"

field_map = open_B_map(path,filename)


# initialize new particle

for it_p in range(0,nb_part):
    it_z = 0
    divX = 0 # 0.05 * np.random.uniform(0,1) 
    divY = 0.05 * np.random.uniform(0,1) 
    
    refp = EtoP(refE)
    dponp = 0
#    px = refp*divX
#    py = refp*divY
#    pz = sqrt(refp**2 - px**2 - py**2) 
#    
    x = 0
    y = 0
    z = 0
#    vx = PtoV(px)
#    vy = PtoV(py)
#    vz = PtoV(pz)    
    beam[0,it_p,:] = [z,x,divX,y,divY,0,dponp]
    
    

    L = 0.15
    [beam[:,it_p,:],it_z] = drift(L,beam[:,it_p,:],refE,it_z,1)

    [beam[:,it_p,:],it_z] = particle_in_B_map(field_map,beam[:,it_p,:],Delta_z,it_z,nb_it_z,it_p,refp,0.002,0.002,spatial_res = 10)

    L = 0.55
    [beam[:,it_p,:],it_z] = drift(L,beam[:,it_p,:],refE,it_z,1)
    
    

print("--- execution time: %s seconds ---" % (time.perf_counter() - start_time))



plt.figure(0)
plt.plot(beam[0:it_z+1,:,0],beam[0:it_z+1,:,1])
plt.figure(1)
plt.plot(beam[0:it_z+1,:,0],beam[0:it_z+1,:,3])


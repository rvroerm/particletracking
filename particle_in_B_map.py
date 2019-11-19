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
from transfer_functions import EOM_particle_in_Bfield, PtoV, EtoP, drift,quad

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

def particle_in_B_map(field_map,beam,itz0,it_p,refp,nb_segments = 10):
    
    
    #List unique values in z
    z_list = field_map['z_m'].unique()
    Delta_z = (max(z_list)-min(z_list))/nb_segments
    
    
    beam_ref = np.empty(shape=[nb_segments + 1,7]) 
    z = beam[itz0,0]
    
    for itz in range(itz0 + 1 , itz0 + nb_segments + 1):    
        
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
        
        
            
        
        beam_ref[itz-itz0,:] = EOM_particle_in_Bfield(Delta_z,[z,0,0,0,0,0,0],refp,Bx,By,Bz)
        #beam_ref[itz-itz0,:] = EOM_particle_in_Bfield_test(Delta_z,[z,0,0,0,0,0,0],refp,Bx,By,Bz)
        
#        # compute rotation coefficients on X axis (case of dipole component)
#        costheta_mat[itz-itz0] = 1/sqrt(1+beam_ref[itz-itz0,2]**2)
#        sintheta_mat[itz-itz0] = beam_ref[itz-itz0,2]/sqrt(1+beam_ref[itz-itz0,2]**2)
#        
        
        # 3) Actual particle
        
            
        if itz == 0:
            print("beam should be initialized before fumction, itz0>0")
        else:
            x = beam[itz-1,1]
            y = beam[itz-1,3]
            
        
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
        
        beam[itz,:] = EOM_particle_in_Bfield(Delta_z,beam[itz-1,:],refp,Bx,By,Bz)
        #beam[itz,:] = EOM_particle_in_Bfield_test(Delta_z,beam[itz-1,:],refp,Bx,By,Bz)
        
        
        # compute difference with reference particle on Xaxis
        beam[itz,1] = beam[itz,1] - beam_ref[itz-itz0,1]
        beam[itz,2] = beam[itz,1] - beam_ref[itz-itz0,2]
        
        #vx = costheta_mat[itz-itz0]*vx - sintheta_mat[itz-itz0]*vz
        #vz = sintheta_mat[itz-itz0]*vx + costheta_mat[itz-itz0]*vz
        
        ## normalize to correct rounding errors
        #v = sqrt(vx**2 + vy**2 + vz**2)
        #vx = vx * v0/v
        #vy = vy * v0/v
        #vz = vz * v0/v
        
        #beam[itz,it_p,:] = [x,y,z,vx,vy,vz]
            
    return [beam,itz]



plt.close('all')




nb_part=10

nb_it_z = 30
Delta_z = 0.01

nb_pts_z = 1000

beam = np.zeros(shape=[nb_pts_z,nb_part,7]) 

refE = 100
refp = EtoP(refE)

path = "D:/beamline/magnets/field maps/test/"
filename = "test_quad_B1T_aperture15_noise0.csv"

field_map = open_B_map(path,filename)


# initialize new particle

for it_p in range(0,nb_part):
    it_z = 0
    divX = 0 #0.05 * np.random.uniform(0,1) 
    divY = 0.005 #* np.random.uniform(0,1) 
    
    print(' ------------  particle ',it_p, '--------------')
    
    dponp = 0
#    px = refp*divX
#    py = refp*divY
#    pz = sqrt(refp**2 - px**2 - py**2) 
#    
    x = 0
    y = 0 #0.05 * np.random.uniform(0,1)
    z = 0
#    vx = PtoV(px)
#    vy = PtoV(py)
#    vz = PtoV(pz)    
    beam[0,it_p,:] = [z,x,divX,y,divY,0,dponp]
    
    L = 0.11
    [beam[:,it_p,:],it_z] = drift(L,beam[:,it_p,:],refE,it_z,1)
    
    #L = 0.0001
    #[beam[:,it_p,:],it_z] = drift(L,beam[:,it_p,:],refE,it_z,1)
    
    use_map = False
    N_segments = 30
    if it_p < nb_part/2:
        [beam[:,it_p,:],it_z] = particle_in_B_map(field_map,beam[:,it_p,:],it_z,it_p,refp,N_segments)
    else:
        L=0.15
        B=-1.
        a=0.015
        [beam[:,it_p,:],it_z] = quad(L,B,a,beam[:,it_p,:],it_z,refE,N_segments)
    
    
    L = 0.5
    [beam[:,it_p,:],it_z] = drift(L,beam[:,it_p,:],refE,it_z,1)
    
    
    print('z = ',beam[it_z,it_p,0])

print("--- execution time: %s seconds ---" % (time.perf_counter() - start_time))


  

plt.figure(0)
plt.plot(beam[0:it_z+1,:,0],beam[0:it_z+1,:,1])
plt.grid(which='major')
plt.figure(1)
plt.plot(beam[0:it_z+1,:,0],beam[0:it_z+1,:,3])
plt.grid(which='major')

plt.figure(3)
plt.scatter(beam[it_z,:,3],beam[it_z,:,4])
plt.xlabel('y [m]')
plt.ylabel('divY [rad]')
plt.grid(which='major')


plt.figure(55)
plt.plot(beam[0:it_z+1,:,0],1/np.vectorize(sqrt)(1+beam[0:it_z+1,:,2]**2+beam[0:it_z+1,:,4]**2))
plt.xlabel('z [m]')
plt.ylabel(' p_z/p_tot')

plt.figure(56)
plt.plot(beam[0:it_z+1,:,0],1/3/10**8*np.vectorize(PtoV)(refp/np.vectorize(sqrt)(1+beam[0:it_z+1,:,2]**2+beam[0:it_z+1,:,4]**2)))
plt.xlabel('z [m]')
plt.ylabel('v_z [1/c]')
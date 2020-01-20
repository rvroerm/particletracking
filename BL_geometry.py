# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:15:05 2020

@author: rvroerm
"""

import numpy as np
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt



def GTR_layout_from_transport(transport_file, nb_pts_z, refE):
    # CCT case; CCT angle in radians
    # aperture_size = vertical aperture
    
    angle = 0
    rot_mat = [[1,0],[0,1]]
    new_vect = [0,0]
    
    beam_coord = np.zeros(shape=[nb_pts_z,2]) # --> initialize with non?
    #beam_coord[0,:] = [0,0]
    
    
    it_z = 1
    
    Brho = 1/300*sqrt(refE**2+2*938*refE)
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       cnt = 1
       
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           
           if data: # string is not empty
               
               if data[0][0:2] == "3.":
                   # drift
                   L = float(data[1].replace(";","") )
                   
                   beam_coord[it_z,:] = beam_coord[it_z-1,:] + np.matmul(rot_mat,[L,0])
                   
                   it_z = it_z + 1
                   
               if data[0][0:2] == "4.":
                   # sector bending
                   
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )
                   rho = Brho/B
                   B_angle = L/rho
                   
                   beam_coord[it_z,:] = beam_coord[it_z-1,:] + np.matmul(rot_mat,[rho*sin(B_angle),rho*(1-cos(B_angle))])
                   
                   
                   angle = angle + B_angle
                   
                   
                   rot_mat = [[cos(angle),-sin(angle)],[sin(angle),cos(angle)]]
                   it_z = it_z + 1
                   
               if data[0][0:2] == "5.": 
                   # quad
                   L = float(data[1].replace(";","") )
                   new_vect[0] = L
                   new_vect[1] = 0
                   beam_coord[it_z,:] = beam_coord[it_z-1,:] + np.matmul(rot_mat,new_vect)
                   it_z = it_z + 1
                   
                   
           line = fp.readline()
           cnt += 1
       


    return beam_coord  









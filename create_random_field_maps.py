# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 10:37:07 2019

@author: rvroerm
"""
import numpy as np
import csv
from math import sqrt

path = "D:/beamline/magnets/field maps/test/"
# create random field maps

step = 2 # mm

x_max = 30
x_min = -30
y_max = 30
y_min = -30
z_max = 150
z_min = 0

noise_level = 0 # 10**-4

with open(path + 'test_quad_B1T_aperture15_noise0.csv','w', newline='') as csvfile:
    
    filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    
    #filewriter.writerow(['Magnet field map generated in python'])
    filewriter.writerow(['x [mm]', 'y [mm]','z [mm]','Bx [T]','By [T]','Bz [T]','dBx/dy [T/m]','dBy/dx [T/m]'])
    
    
    
    
    for it_Z in range(0,int((z_max-z_min)/step) + 1):
        
        z = z_min + it_Z*step
        for it_Y in range(0,int((y_max-y_min)/step) + 1):
            y = y_min + it_Y*step
            for it_X in range(0,int((x_max-x_min)/step) + 1):
                
                x = x_min + it_X*step
                
                dBx_dy = - 1/0.015
                dBy_dx = - 1/0.015
                
                Bx = dBx_dy * y/1000 + np.random.normal(0,noise_level)
                By = dBy_dx * x/1000 + np.random.normal(0,noise_level)
                Bz = np.random.normal(0,noise_level)
                
                
                
                filewriter.writerow([x,y,z,Bx,By,Bz,dBx_dy,dBy_dx])
            
               
                

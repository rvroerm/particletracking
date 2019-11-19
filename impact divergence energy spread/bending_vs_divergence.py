# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:40:08 2019

@author: HILAM
"""
from math import sqrt,sin,cos
import numpy as np
import matplotlib.pyplot as plt


alpha = 0.05
rot_mat1 = [[cos(alpha) , -sin(alpha)],[sin(alpha) , cos(alpha)]]
rot_mat2 = [[cos(alpha) , sin(alpha)],[-sin(alpha) , cos(alpha)]]

x = np.arange(0,0.5,0.001)
nb_pts = len(x)
x1 = np.zeros(nb_pts)
y1 = np.zeros(nb_pts)
x2 = np.zeros(nb_pts)
y2 = np.zeros(nb_pts)

B=1.5


E=200
Brho  = 1/300*sqrt(E**2+2*938*E)
R=Brho/B

y = np.sqrt(R**2 - (x-R)**2)


# rotation for beam spread
for i in range(0,nb_pts):
    [x1[i] , y1[i]] = np.matmul(rot_mat1,[x[i],y[i]])
    [x2[i] , y2[i]] = np.matmul(rot_mat2,[x[i],y[i]])

plt.figure(1)
plt.plot(x1,y1,color='blue')
plt.plot(x2,y2,color='blue')

#####

# other energy
x1 = np.zeros(nb_pts)
y1 = np.zeros(nb_pts)
x2 = np.zeros(nb_pts)
y2 = np.zeros(nb_pts)

E=180
Brho  = 1/300*sqrt(E**2+2*938*E)
R=Brho/B

y = np.sqrt(R**2 - (x-R)**2)

# rotation for beam spread
for i in range(0,nb_pts):
    [x1[i] , y1[i]] = np.matmul(rot_mat1,[x[i],y[i]])
    [x2[i] , y2[i]] = np.matmul(rot_mat2,[x[i],y[i]])

plt.plot(x1,y1,color='red')
plt.plot(x2,y2,color='red')
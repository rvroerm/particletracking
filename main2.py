# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:18:03 2020

@author: rvroerm
"""

from BL_classes import create_BL_from_Transport, BL_plot, Particle, Beam
from math import pi
import matplotlib.pyplot as plt
import time


plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = pi/6)


[fig, ax_X, ax_Y] = BL_plot(my_beamline)

my_beam = Beam(nb_part=10, refE = 160, DeltaE=5, E_dist='uniform',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='cst')




for particle in my_beam.particle_list :
    
    t0 = time.time()
    particle.particle_through_BL(my_beamline)
    print("exec time tot  %.2E \n "%(time.time() - t0))
    
    ax_X.plot(particle.z_df.values, particle.p_df[['x [m]']].values)
    ax_Y.plot(particle.z_df.values, particle.p_df[['y [m]']].values)
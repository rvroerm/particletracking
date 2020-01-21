# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:18:03 2020

@author: rvroerm
"""

from BL_classes import create_BL_from_Transport, BL_plot, Particle, Beam
from math import pi
import matplotlib.pyplot as plt



plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = pi/6)


#BL_plot(my_beamline, orientation='ZY')

my_beam = Beam(nb_part=2, refE = 160, DeltaE=0, E_dist='cst',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')


plt.figure('test')

for particle in my_beam.beam_df['Particle objects'].values :
    
    particle.particle_through_BL(my_beamline)
    
    plt.plot(particle.z_df.values, particle.p_df[['x [m]']].values)
    
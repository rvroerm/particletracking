# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:18:03 2020

@author: rvroerm
"""

from BL_classes import create_BL_from_Transport, BL_plot, Particle, Beam
from math import pi
import matplotlib.pyplot as plt
import time
import numpy as np

plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = pi/6)



#########################################################################
# plot beam through BL

[fig, ax_X, ax_Y] = BL_plot(my_beamline)

my_beam = Beam(nb_part=1000, refE = 160, DeltaE=1.5, E_dist='uniform',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='cst')



t0 = time.time()


for particle in my_beam.particle_list :
    
    
    particle.particle_through_BL(my_beamline)
    
    
    ax_X.plot(particle.z[0:particle.it], particle.X[0:particle.it,0])
    ax_Y.plot(particle.z[0:particle.it], particle.X[0:particle.it,2])
    
    
    
print("exec time tot  %.2E \n "%(time.time() - t0))



###############################################################################\

# transverse plot at ISO

plt.figure('transverse profile')


index = my_beamline.get_element_index("ISO")
z_ISO = my_beamline.BL_df.loc[index,'z [m]']

p_index = np.where( abs(my_beam.particle_list[0].z[0:particle.it] - z_ISO) < 10**-6 )
if len(p_index[0]) == 1 :
    p_index = np.asscalar(p_index[0])
else:
    raise Exception('ISO not found')

X = my_beam.get_beam_x(p_index)

efficiency = len(X[~np.isnan(X)]) / len(X)
print("efficiency = %.2f %% "%(efficiency*100))

nb_bins = 20
plt.hist(X,nb_bins,range=(-0.03,0.03))

###############################################################################
# tune BL element

#my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)
#
#[fig, ax_X, ax_Y] = BL_plot(my_beamline)
#
#
#element_to_tune = "Q4X"
#plt.suptitle("Tuning of %s"%element_to_tune)
#
#index = my_beamline.get_element_index(element_to_tune)
#ref_field = my_beamline.BL_df.loc[index,'BL object'].Bfield
#tune_range = [0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.1]
#
#for B_factor in tune_range:
#    my_proton = Particle(z=0, x=0, y=0, divX=0.025, divY=0.025, p=570.75, refp=570.75, max_it=1000)
#    
#    # change field
#    my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field  * B_factor
#    
#    
#    my_proton.particle_through_BL(my_beamline)
#    
#    ax_X.plot(my_proton.z[0:my_proton.it], my_proton.X[0:my_proton.it,0])
#    ax_Y.plot(my_proton.z[0:my_proton.it], my_proton.X[0:my_proton.it,2])
#    
#    
#
## reset field    
#my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field    
#    
#ax_X.legend([(element_to_tune + " = " + str(round(x*ref_field,2))) for x in tune_range], loc="upper right")
#ax_Y.legend([(element_to_tune + " = " + str(round(x*ref_field,2))) for x in tune_range], loc="upper right") 




    



   
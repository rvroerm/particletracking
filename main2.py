# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 13:18:03 2020

@author: rvroerm
"""

from BL_classes import create_BL_from_Transport, Particle, Beam
from math import pi
import matplotlib.pyplot as plt
import seaborn as sns
import time
import numpy as np
from scipy.stats import norm
from transfer_functions import EtoP
from plot_functions import plot_beam_through_BL, BL_plot

plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = pi/6)



#########################################################################
# plot beam through BL

refE = 160

my_beam = Beam(nb_part=1000, refE = refE, DeltaE=5, E_dist='cst',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')

t0 = time.time()

[fig, ax_X, ax_Y] = plot_beam_through_BL(my_beam = my_beam, my_beamline = my_beamline)
  
    
print("exec time tot  %.2E \n "%(time.time() - t0))



###############################################################################\

# transverse plot at ISO


index = my_beamline.get_element_index("ISO")
z_ISO = my_beamline.BL_df.loc[index,'z [m]']

p_ISO_index = my_beam.particle_list[0].get_z_index(z_ISO)


#X = my_beam.get_beam_x(p_ISO_index)

X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
Y = my_beam.get_beam_param(param='y', row_nb=p_ISO_index)

efficiency = len(X[~np.isnan(X)]) / len(X)
print("efficiency = %.2f %% "%(efficiency*100))


fig = plt.figure('transverse profile',figsize=(12, 6))
    
# increase space between plots
fig.subplots_adjust(wspace=0.4)

ax_X = fig.add_subplot(1,2,1)
ax_Y = fig.add_subplot(1,2,2)

sns.distplot(np.clip(X[~np.isnan(X)],-0.1,0.1)*1000, kde=False, fit=norm, ax=ax_X)
(mu, sigma) = norm.fit(X[~np.isnan(X)])
ax_X.legend(["sigma = %.2f mm"%(sigma*1000)], loc='upper right')
ax_X.set_xlabel('x [mm]')

sns.distplot(np.clip(Y[~np.isnan(Y)],-0.1,0.1)*1000, kde=False, fit=norm, ax=ax_Y)
(mu, sigma) = norm.fit(Y[~np.isnan(Y)])
ax_Y.legend(["sigma = %.2f mm"%(sigma*1000)], loc='upper right')
ax_Y.set_xlabel('y [mm]')



###############################################################################\
# energy plot

fig = plt.figure("energy plot")

E_source = my_beam.get_beam_param(param='E', row_nb=0)
E_ISO = my_beam.get_beam_param(param='E', row_nb=p_ISO_index)
plt.hist(E_source[~np.isnan(E_source)], bins=20, alpha=0.3, range=(refE-5,refE+5), label="E source")
plt.hist(E_ISO[~np.isnan(E_ISO)], bins=20, alpha=0.3, range=(refE-5,refE+5), label="E ISO")
plt.title('Energy transmission')
plt.legend(loc='upper right')
plt.xlabel('E [MeV]')
plt.ylabel('nb of protons (arb. units)')



###############################################################################
# tune BL element

#my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)
#
#[fig, ax_X, ax_Y] = BL_plot(my_beamline)
#
#
#element_to_tune = "Q5Y"
#plt.suptitle("Tuning of %s"%element_to_tune)
#
#index = my_beamline.get_element_index(element_to_tune)
#ref_field = my_beamline.BL_df.loc[index,'BL object'].Bfield
#tune_range = [0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.1]
#tune_range = [0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03]
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
#ax_X.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right")
#ax_Y.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right") 




    



   
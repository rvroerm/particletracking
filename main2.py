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
from plot_functions import BL_geometry, plot_beam_through_BL, BL_plot_for_traces


t0 = time.time()

plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)



refE = 160


# plot beamline
BL_geometry(my_beamline, refp=EtoP(refE))

#########################################################################
# plot beam through BL


my_beam = Beam(nb_part=100, refE = refE, DeltaE=1.5, E_dist='cst',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')



[fig, ax_X, ax_Y] = plot_beam_through_BL(my_beam = my_beam, my_beamline = my_beamline)
  
    




###############################################################################\

# transverse plot at ISO


index = my_beamline.get_element_index("ISO")
z_ISO = my_beamline.BL_df.loc[index,'z [m]']

p_ISO_index = my_beam.particle_list[0].get_z_index(z_ISO)

print(index)
print(p_ISO_index)
print("z = ",my_beam.get_beam_param(param='z', row_nb=p_ISO_index)[0])

#X = my_beam.get_beam_x(p_ISO_index)

X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
Y = my_beam.get_beam_param(param='y', row_nb=p_ISO_index)

efficiency = len(X[~np.isnan(X)]) / len(X)
print("efficiency = %.2f %% "%(efficiency*100))


fig = plt.figure('transverse profile',figsize=(9, 9))

sns.distplot(np.clip(X[~np.isnan(X)],-0.1,0.1)*1000, kde=False, fit=norm, fit_kws={"color":"blue"})
(muX, sigmaX) = norm.fit(X[~np.isnan(X)])

sns.distplot(np.clip(Y[~np.isnan(Y)],-0.1,0.1)*1000, kde=False, fit=norm, fit_kws={"color":"red"})
(muY, sigmaY) = norm.fit(Y[~np.isnan(Y)])

plt.legend(["sigmaX = %.2f mm"%(sigmaX*1000),"sigmaY = %.2f mm"%(sigmaY*1000)], loc='upper right')

plt.xlabel('x/y [mm]')
plt.ylabel('counts (arb units)')



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

my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)

[fig, ax_X, ax_Y] = BL_plot_for_traces(my_beamline, title='variation study 1 particle')

element_to_tune = "Q5Y"
plt.suptitle("Tuning of %s"%element_to_tune)

index = my_beamline.get_element_index(element_to_tune)
ref_field = my_beamline.BL_df.loc[index,'BL object'].Bfield
tune_range = [0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.1]
#tune_range = [0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03]

for B_factor in tune_range:
    
    ##############################
    # 1. case of 1 particle at ref E
    
    my_proton = Particle(z=0, x=0, y=0, divX=0.025, divY=0.025, p=570.75, refp=570.75, max_it=1000)
    
    # change field
    my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field  * B_factor
    
    
    my_proton.particle_through_BL(my_beamline)
    
    ax_X.plot(my_proton.z[0:my_proton.it], my_proton.X[0:my_proton.it,0])
    ax_Y.plot(my_proton.z[0:my_proton.it], my_proton.X[0:my_proton.it,2])
    
    
    ##############################
    # 2. case of a full beam
    
    my_beam = Beam(nb_part=100, refE = refE, DeltaE=1.5, E_dist='cst',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')
    
    my_beam.beam_through_BL(my_beamline)
    sigma_X = my_beam.size_X(row_nb = p_ISO_index)
    sigma_Y = my_beam.size_Y(row_nb = p_ISO_index)
    print('B_factor = %0.2f : sigma_X = %0.2f and sigma_Y = %0.2f '%(B_factor, sigma_X*1000, sigma_Y*1000))
    
    

# reset field    
my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field    
    
ax_X.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right")
ax_Y.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right") 




    
print("exec time tot  %.2E \n "%(time.time() - t0))


   
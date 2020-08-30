# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 10:53:58 2020

@author: rvroerm
"""

# test impact of source offset versus the beamline



from BL_classes import create_BL_from_Transport, Particle, Beam
from math import pi
import matplotlib.pyplot as plt
import seaborn as sns
import time
import numpy as np
from scipy.stats import norm
from transfer_functions import EtoP
from plot_functions import BL_geometry, plot_beam_through_BL, BL_plot_for_traces

import pandas as pd

t0 = time.time()

plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)



refE = 160
DeltaE = 1.6

# plot beamline
BL_geometry(my_beamline, refp=EtoP(refE))



index = my_beamline.get_element_index("ISO")
z_ISO = my_beamline.BL_df.loc[index,'z [m]']



#########################################################################
# plot beam through BL

max_offset = 0.002
offset_list = np.arange(-max_offset, max_offset + max_offset/10, max_offset/10)

results_df = pd.DataFrame(columns = ['offset [mm]','BL eff [%]','size X [mm]','size Y [mm]'])

for offset in offset_list :

    my_beam = Beam(nb_part=1000, refE = refE, DeltaE=DeltaE, E_dist='uniform2',  \
                            DeltaX = 10**-5, DeltaY = 10**-5, size_dist='uniform', \
                            DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform', \
                            OffsetX = offset, OffsetY = offset)
    
    
        
    my_beam.beam_through_BL(my_beamline)
    
    p_ISO_index = my_beam.particle_list[0].get_z_index(z_ISO)
    
    
    X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
    efficiency = len(X[~np.isnan(X)]) / len(X)
    
    sigma_X = my_beam.size_X(row_nb = p_ISO_index)
    sigma_Y = my_beam.size_Y(row_nb = p_ISO_index)
  
    #print("offset = %.2f mm : efficiency = %.2f %%, sigmaX= %.2f mm and sigmaY= %.2f mm"%(offset*1000, efficiency*100, sigma_X*1000, sigma_Y*1000))
    
    new_line = pd.DataFrame.from_dict({'offset [mm]' : [offset*1000], 
                                        'BL eff [%]' : [efficiency*100], 
                                        'size X [mm]' : [sigma_X*1000], 
                                        'size Y [mm]': [sigma_Y*1000]})
    results_df = results_df.append(new_line, ignore_index=True)
    
    
print(results_df.round({'offset [mm]':2, 'BL eff [%]':0, 'size X [mm]':2, 'size Y [mm]':2}))
    

    


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

import pandas as pd

t0 = time.time()

plt.close('all')

# parameters
input_file = "C:/TRANS/for001.dat"
my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)



refE = 160
DeltaE = 2

# plot beamline
BL_geometry(my_beamline, refp=EtoP(refE))






#########################################################################
# plot beam through BL


my_beam = Beam(nb_part=1000, refE = refE, DeltaE=DeltaE, E_dist='uniform2',  \
                        DeltaX = 10**-5, DeltaY = 10**-5, size_dist='normal', \
                        DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform', \
                        OffsetX = 0., OffsetY=0)


# my_beam = Beam(nb_part=1, refE = refE, DeltaE=0, E_dist='cst',  \
#                         DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
#                         DeltaDivX = 0.025, DeltaDivY = 0.025, div_dist='cst')

    
[fig, ax_X, ax_Y] = plot_beam_through_BL(my_beam = my_beam, my_beamline = my_beamline)
  
    




###############################################################################\

# transverse plot at ISO


index = my_beamline.get_element_index("ISO")
z_ISO = my_beamline.BL_df.loc[index,'z [m]']

p_ISO_index = my_beam.particle_list[0].get_z_index(z_ISO)

print("z ISO = ",my_beam.get_beam_param(param='z', row_nb=p_ISO_index)[0])

#X = my_beam.get_beam_x(p_ISO_index)

X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
Y = my_beam.get_beam_param(param='y', row_nb=p_ISO_index)

efficiency = len(X[~np.isnan(X)]) / len(X)
print("efficiency = %.2f %% "%(efficiency*100))


fig = plt.figure('transverse profile', figsize=(9, 9))
ax = fig.add_subplot(111)

clip_dist = 0.2 # show particles within [-clip_dist,clip_dist] only

sns.distplot(np.clip(X[~np.isnan(X)], -clip_dist, clip_dist)*1000, kde=False, fit=norm, fit_kws={"color":"blue"}, ax=ax)
sns.distplot(np.clip(Y[~np.isnan(Y)], -clip_dist, clip_dist)*1000, kde=False, fit=norm, fit_kws={"color":"red"}, ax=ax)


posX = my_beam.pos_X(row_nb = p_ISO_index, clip_dist=clip_dist)
posY = my_beam.pos_Y(row_nb = p_ISO_index, clip_dist=clip_dist)
print("posX = %.2f mm, posY = %.2f mm"%(posX*1000 , posY*1000))

sigmaX = my_beam.size_X(row_nb = p_ISO_index, clip_dist=clip_dist)
sigmaY = my_beam.size_Y(row_nb = p_ISO_index, clip_dist=clip_dist)
print("sigmaX = %.2f mm, sigmaY = %.2f mm"%(sigmaX*1000 , sigmaY*1000))
ax.legend(["sigmaX = %.2f mm"%(sigmaX*1000),"sigmaY = %.2f mm"%(sigmaY*1000)], loc='upper right')


ax.set_xlabel('x/y [mm]')
ax.set_ylabel('counts (arb units)')


# plot X-Y positions at isocenter with color legend refering to the particle energy
E = my_beam.get_beam_param(param='E', row_nb=0)

fig = plt.figure('positions vs E')
sns.scatterplot(x=X, y=Y, hue=E, hue_norm=(refE-abs(DeltaE), refE+abs(DeltaE)))
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
plt.legend(title = 'E [MeV]')


###############################################################################\
# energy plot

fig = plt.figure("energy plot")

E_source = my_beam.get_beam_param(param='E', row_nb=0)
E_ISO = my_beam.get_beam_param(param='E', row_nb=p_ISO_index)
plt.hist(E_source[~np.isnan(E_source)], bins=20, alpha=0.3, range=(refE-abs(DeltaE), refE+abs(DeltaE)), label="E source")
plt.hist(E_ISO[~np.isnan(E_ISO)], bins=20, alpha=0.3, range=(refE-abs(DeltaE), refE+abs(DeltaE)), label="E ISO")
plt.title('Energy transmission')
plt.legend(loc='upper right')
plt.xlabel('E [MeV]')
plt.ylabel('nb of protons (arb. units)')



###############################################################################
# tune BL element

variation_study = False

varaiation_df = pd.DataFrame(columns = ['magnet','factor','size X','size Y','BL eff'])

if variation_study:
    
    print('\n')
    
    my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)
    
    elements_to_tune = ["Q1C","Q2C","Q3C","Q4X","Q5Y","Q1G","Q2G","Q3G","Q4G"]
    elements_to_tune = ["Q4X","Q5Y"]
    elements_to_tune = ["TUN1","TUN2"]
    
    #tune_range = [0.95, 0.975, 0.99, 1, 1.01, 1.025, 1.05]
    tune_range = [0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2]
    
    
    for element_to_tune in elements_to_tune:
        print("Element: %s"%element_to_tune)
        
        # get element from beamline
        index = my_beamline.get_element_index(element_to_tune)
        ref_field = my_beamline.BL_df.loc[index,'BL object'].Bfield
        
        [fig, ax_X, ax_Y] = BL_plot_for_traces(my_beamline, title='variation study 1 particle, element %s'%element_to_tune)
        plt.suptitle("Tuning of %s"%element_to_tune)
        
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
            
            my_beam = Beam(nb_part=100, refE = refE, DeltaE=1.5, E_dist='uniform2',  \
                               DeltaX = 10**-5, DeltaY = 10**-5, size_dist='uniform', \
                               DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')
            
            my_beam.beam_through_BL(my_beamline)
            
            p_ISO_index = my_beam.particle_list[0].get_z_index(z_ISO)
            sigma_X = my_beam.size_X(row_nb = p_ISO_index)
            sigma_Y = my_beam.size_Y(row_nb = p_ISO_index)
            
            # efficiency
            X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
            
            efficiency = len(X[~np.isnan(X)]) / len(X)
            print('Field = %0.2f'%(ref_field * B_factor))
            print('B_factor = %0.2f : sigma_X = %0.2f , sigma_Y = %0.2f , efficiency = %.1f %% '%(B_factor, sigma_X*1000, sigma_Y*1000, efficiency*100))
            

            
            new_line = pd.DataFrame.from_dict({'magnet' : [element_to_tune], 
                                               'factor' : [B_factor], 
                                               'size X' : [sigma_X*1000], 
                                               'size Y': [sigma_Y*1000],
                                               'BL eff': [efficiency]})
            varaiation_df = varaiation_df.append(new_line, ignore_index=True)
            
            
    
        # reset field    
        my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field    
        
        ax_X.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right")
        ax_Y.legend([(element_to_tune + " = " + str(round(x*ref_field,3))) for x in tune_range], loc="upper right")
    
        print('\n')



###############################################################################
# correcting magnet variation

variation_study = False
clip_dist = 0.03

if variation_study:
    
    
    
    my_beamline = create_BL_from_Transport(input_file, CCT_angle = 0)
    
    #elements_to_tune = ["Q1C","Q2C","Q3C","Q4X","Q5Y","Q1G","Q2G","Q3G","Q4G"]
    
    elements_to_correct = ["Q4X"]
    #elements_to_correct = ["Q4X","Q5Y"]
    error_range = [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.5, 2]
    error_range = [0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2]
    
    elements_to_tune = ["QN1", "QN3"]
    elements_to_tune = ["Q5Y"]
    
    #tune_range = [0.95, 0.975, 0.99, 1, 1.01, 1.025, 1.05]
    #tune_range = [0.5, 0.75, 0.9, 0.95, 1, 1.05, 1.1, 1.25, 1.5]
    tune_range = [0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2]
    #tune_range = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.5,2]
    
    for element_to_correct in elements_to_correct:
        
        print('\n')
        print('****************************************')
        print("Element to correct: %s"%element_to_correct)
        
        # get element from beamline
        index1 = my_beamline.get_element_index(element_to_correct)
        
        ref_field1 = my_beamline.BL_df.loc[index1,'BL object'].Bfield
        
        for error_on_magnet in error_range:
        
            my_beamline.BL_df.loc[index1,'BL object'].Bfield = ref_field1  * error_on_magnet
            print('\n %s = %0.2f'%(element_to_correct, ref_field1  * error_on_magnet))
            
            for element_to_tune in elements_to_tune:
                print("Element: %s"%element_to_tune)
                
                # get element from beamline
                index = my_beamline.get_element_index(element_to_tune)
                ref_field = my_beamline.BL_df.loc[index,'BL object'].Bfield
                
                
                for B_factor in tune_range:
                    
                    # change field
                    my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field  * B_factor
                    
                    
                    # case of a full beam
                    
                    my_beam = Beam(nb_part=100, refE = refE, DeltaE=3, E_dist='uniform2',  \
                                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='uniform', \
                                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform')
                    
                    my_beam.beam_through_BL(my_beamline)
                    sigma_X = my_beam.size_X(row_nb = p_ISO_index, clip_dist=clip_dist)
                    sigma_Y = my_beam.size_Y(row_nb = p_ISO_index, clip_dist=clip_dist)
                    
                    # efficiency
                    X = my_beam.get_beam_param(param='x', row_nb=p_ISO_index)
                    
                    efficiency = len(X[~np.isnan(X)]) / len(X)
                    print('B_factor = %0.2f : sigma_X = %0.2f , sigma_Y = %0.2f , efficiency = %.1f %% '%(B_factor, sigma_X*1000, sigma_Y*1000, efficiency*100))
                
                
                
        
                # reset field    
                my_beamline.BL_df.loc[index,'BL object'].Bfield = ref_field    
            
        
        # reset field    
        my_beamline.BL_df.loc[index1,'BL object'].Bfield = ref_field1      
        
    print('\n')


    
print("exec time tot  %.2E \n "%(time.time() - t0))


"""
Created on Tue Jan 21 13:18:03 2020
@author: rvroerm
"""

from BL_classes import create_BL_from_Transport, Particle, Beam, Beamline, BL_Element, ScanMagnet
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
DeltaE = 5



#########################################################################
# plot beam through BL


my_beam = Beam(nb_part=10, refE = refE, DeltaE=DeltaE, E_dist='uniform2',  \
                        DeltaX = 10**-5, DeltaY = 10**-5, size_dist='uniform', \
                        DeltaDivX = 0.0005, DeltaDivY = 0.0005, div_dist='uniform')
    
    
my_beamline = Beamline()

drift = BL_Element(length = 0.2)

SM = ScanMagnet(By = 0.3, length=0.2)


my_beamline.add_element(drift)
my_beamline.add_element(SM)
my_beamline.add_element(BL_Element(length = 2))
my_beamline.add_element(drift)

[fig, ax_X, ax_Y] = plot_beam_through_BL(my_beam = my_beam, my_beamline = my_beamline)

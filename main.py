# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 12:04:20 2020

@author: rvroerm
"""

from run_from_transport import run_from_transport


# parameters
input_file = "C:/TRANS/for001.dat"

nb_part=1000
refE = 160
old_refE = 160
DeltaE=5
E_dist='cst'

DeltaDivX = 0.05
DeltaDivY = 0.05
div_dist='uniform'
gap = 0.03

paraxial_correction = False
dpzTolerance = 10**-4



[sigmaX, sigmaY, eff_ESS_dEonE_1pc, eff_GTR_dEonE_1pc] = \
run_from_transport(input_file = input_file, nb_part=nb_part, \
                       N_segments = 10, kill_lost_particles = True, \
                       refE = refE, old_refE = old_refE, DeltaE=DeltaE, E_dist=E_dist,  \
                       DeltaDivX = DeltaDivX, DeltaDivY = DeltaDivY, div_dist=div_dist, \
                       gap = gap, paraxial_correction = False, \
                       plot_results = True, output_results=True)
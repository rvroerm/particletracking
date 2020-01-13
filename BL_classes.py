# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""

def sec_bend(L,B,n,gap,beam,refE,it_z,N_segments = 10,kill_lost_particles=True,rest_mass=938,gap_X=1, \
          paraxial_correction = False, dpz_tolerance=10**-6):


class Dipole:
    
    def __init__(self, length, field, order, gap, pole_face1, pole_face2):
        self.L = length
        self.B = field
        self.n = order
        
    def setField(self, B):
        self.B = B
        
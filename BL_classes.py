# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""
import numpy as np
from scipy.stats import norm

class BeamParticle:
    
    def __init__(self, x=0, y=0, divX=0, divY=0, energy_MeV=0, particle_type='proton'):
        self.x = x
        self.y = y
        self.divX = divX
        self.divY = divY
        self.E = energy_MeV
        if particle_type=='proton' : self.mass = 938
        else: raise Exception("Particle type is not known")
    


class Beam:
    
    def __init__(self):
        self.X = []
        self.Y = []
        self.divX = []
        self.divY = []
        self.E = []
        
    def add_particle(self, particle : BeamParticle):
        self.X = self.X.append(particle.x)
        self.Y = self.Y.append(particle.y)
        self.divX = self.divX.append(particle.divX)
        self.divY = self.divY.append(particle.divY)
        self.E = self.E.append(particle.E)
        
    def size_X(self):
        (muX, sigmaX) = norm.fit(self.X[~np.isnan(self.X)])
        return sigmaX
    
    def size_Y(self):
        (muY, sigmaY) = norm.fit(self.Y[~np.isnan(self.Y)])
        return sigmaY
    


    
class BL_element:
    """
    Define any element in the beamline.  Default = drift
    """

    def __init__(self, name="", length=0):
        self.name = name
        self.length = length
        self.element_type = "drift" # default is a drift
        self.N_segments = 1
        
    def dipole(self, B, n, gap, pole_face1=0, pole_face2=0, nb_pts=10):
        self.element_type = "dipole"
        self.Bfield = B
        self.order = n
        self.gap = gap
        self.N_segments = nb_pts
        
    def quadrupole(self, B, a, nb_pts=10):
        self.element_type = "quad"
        self.Bfield = B
        self.aperture = a
        self.N_segments = nb_pts
        
    def sextupole(self, B, a, nb_pts=10):
        self.element_type = "sext"
        self.Bfield = B
        self.aperture = a
        self.N_segments = nb_pts
        
    def solenoid(self, B, a, nb_pts=20):
        self.element_type = "solenoid"
        self.Bfield = B
        self.aperture = a
        self.N_segments = nb_pts
        
    def slit(self, gap, orientation='X', mat ='tantalum', offset = 0):
        self.element_type = "slit"
        self.gap = gap
        self.orientation = orientation
        self.material = mat
        self.offset = offset
        
        
        
    



    
class Beamline:
    
    def __init__(self):
        self.z = [0]
        self.elements_list = {}
    
    def add_element(self, element : BL_element):
        
        slice_element = element # divide element in many slices
        slice_element.length = element.length/element.N_segments
        
        for i in np.arange(0, element.N_segments):
            self.elements_list.append(slice_element)
            self.z.append(self.z[-1] + slice_element.length)
        
    
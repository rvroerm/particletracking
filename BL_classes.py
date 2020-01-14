# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from transfer_functions import PtoGamma, PtoBrho, EtoP, PtoE, quad_matrix
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt


    
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
        
        
        self.BL_df = pd.DataFrame([[ 0, "souce", "", BL_element()]], \
                                  columns = ['z [m]', 'Element name', 'Element type', 'BL object'])
        
        self.z = 0
        self.elements_list = []
        
        
        
    def get_z(self, row_nb=-1): return self.BL_df[['z [m]']].iloc[row_nb,0]
    
    def add_element(self, element: BL_element()):
        
        # divide element in many slices
        slice_element = element 
        slice_element.length = element.length/element.N_segments
        
        for i in np.arange(0, element.N_segments):
            
            
            new_row = pd.DataFrame(data = [[self.get_z() + slice_element.length, \
                                            slice_element.name, \
                                            slice_element.element_type, \
                                            slice_element ]], \
                                     columns = self.BL_df.columns)
            
            
            
            self.BL_df = self.BL_df.append(new_row, ignore_index=True)
            
        



class Particle:
    
    def __init__(self, z=0, x=0, y=0, divX=0, divY=0, dL=0, p=0, refp=1, particle_type='proton'):
        
        
        if particle_type=='proton' : self.rest_mass = 938
        else: raise Exception("Particle type is not known")
        
        if ~refp> 0:
            raise Exception("Reference energy is not a positive number")
        
        
        
        # build dataframe with z coordinate
        self.z_df = pd.DataFrame([[z]] , columns = ['z [m]'])
        
        # build dataframe with particle properties
        self.p_df = pd.DataFrame([[x, divX, y, divY, dL, (p-refp)/refp]], \
                            columns = ['x [m]','div x [rad]','y [m]','div y [rad]','dL [m]','dp/p'])
        
        
        
        
        self.refp = refp
        self.particle_type = particle_type
        
    
    
    def get_z(self, row_nb=-1): return self.z_df[['z [m]']].iloc[row_nb,0]
    def get_x(self, row_nb=-1): return self.p_df[['x [m]']].iloc[row_nb,0]
    def get_dponp(self, row_nb=-1): return self.p_df[['dp/p']].iloc[row_nb,0]
    
    
    
    def get_vect(self, row_nb=-1): 
        # get beam properties at row=row_nb, default=last row
        return self.p_df[['x [m]','div x [rad]','y [m]','div y [rad]','dL [m]','dp/p']].iloc[row_nb].to_numpy()
    
    def get_p(self, row_nb=-1): return self.get_dponp(row_nb)*self.refp + self.refp
    
    def get_E(self, row_nb=-1): return PtoE(self.get_p(row_nb))
    
    
        
    
    def particle_through_BLelement(self, element : BL_element):
        
        if np.isnan(self.get_dponp()) :
            # add rows with NAN
            self.p_df = self.p_df.append(pd.Series(), ignore_index=True)
            self.z_df = self.z_df.append(pd.Series(), ignore_index=True)
            return
        
        if element.element_type == "quad" :
            
            B = element.Bfield
            L = element.length
            a = element.aperture
            
            quad_mat = quad_matrix(L, B, a, self.get_p(), rest_mass=self.rest_mass)
            
            
            # calculate new particle properties and put them in a new raw
            new_row = pd.DataFrame(data = [np.matmul(quad_mat, self.get_vect())], \
                                     columns = self.p_df.columns)
            self.p_df = self.p_df.append(new_row, ignore_index=True)
            
            # update z
            self.z_df = self.z_df.append(pd.DataFrame(data = [self.get_z() + L], columns = ['z [m]']), ignore_index=True)
            


    


class Beam:
    
    def __init__(self):
        self.X = []
        self.Y = []
        self.divX = []
        self.divY = []
        self.E = []
        
    def add_particle(self, particle : Particle()):
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
    


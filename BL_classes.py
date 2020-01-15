# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""
import numpy as np
import pandas as pd
from scipy.stats import norm
from transfer_functions import PtoGamma, PtoBrho, EtoP, PtoE, drift_matrix, quad_matrix, sec_bend_matrix, pole_face_mat
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt
import warnings

import matplotlib.pyplot as plt

    
class BL_Element:
    """
    Define any element in the beamline.  Default = drift
    """

    def __init__(self, name="", length=0):
        self.name = name
        self.length = length
        self.element_type = "drift" # default is a drift
        self.N_segments = 1
        
        
        
    def dipole(self, B, n, gap, pole_face1=0, pole_face2=0, k1 = 0.5, k2 = 0, nb_pts=10):
        self.element_type = "dipole"
        self.Bfield = B
        self.order = n
        self.gap = gap
        self.k1 = k1
        self.k2 = k2
        self.angle1 = pole_face1
        self.angle2 = pole_face2
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
        
    def slit(self, gap, orientation='X', mat ='tantalum', offset = 0, nb_pts=20):
        self.element_type = "slit"
        self.gap = gap
        self.orientation = orientation
        self.material = mat
        self.offset = offset
        self.N_segments = nb_pts
        
        

class Beamline:
    
    def __init__(self):
        
        
        self.BL_df = pd.DataFrame([[ 0, "souce", "", BL_Element()]], \
                                  columns = ['z [m]', 'Element name', 'Element type', 'BL object'])
        
        
        self.elements_list = []
        
        
        
    def get_z(self, row_nb=-1): return self.BL_df[['z [m]']].iloc[row_nb,0]
    
    def add_element(self, element: BL_Element()):
        
        # divide element in many slices
        slice_element = element 
        slice_element.length = element.length/element.N_segments
        
        for i in np.arange(0, element.N_segments):
            
            if i!=0:
                # not first sub-element
                if element.element_type == 'dipole':
                    slice_element.angle1 = np.nan
            
            if i != element.N_segments-1:
                # not last sub-element
                if element.element_type == 'dipole':
                    slice_element.angle2 = np.nan
            
            
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
    
    def plot_x(self): # example, to put in separate function later
        plt.plot(self.z_df, self.p_df[['x [m]']])
    
    def get_vect(self, row_nb=-1): 
        # get beam properties at row=row_nb, default=last row
        return self.p_df[['x [m]','div x [rad]','y [m]','div y [rad]','dL [m]','dp/p']].iloc[row_nb].to_numpy()
    
    def get_p(self, row_nb=-1): return self.get_dponp(row_nb)*self.refp + self.refp
    
    def get_E(self, row_nb=-1): return PtoE(self.get_p(row_nb))
    
    
        
    
    def particle_through_BLelement(self, element : BL_Element):
        
        if np.isnan(self.get_dponp()) :
            # add rows with NAN
            self.p_df = self.p_df.append(pd.Series(), ignore_index=True)
            self.z_df = self.z_df.append(pd.Series(), ignore_index=True)
            return
        
        
        if element.element_type == "drift" :
            L = element.length
            drift_mat = drift_matrix(L, self.get_p(), rest_mass=self.rest_mass)
            
            # calculate new particle properties and put them in a new raw
            new_row = pd.DataFrame(data = [np.matmul(drift_mat, self.get_vect())], \
                                     columns = self.p_df.columns)
            self.p_df = self.p_df.append(new_row, ignore_index=True)
            
            # update z
            self.z_df = self.z_df.append(pd.DataFrame(data = [self.get_z() + L], columns = ['z [m]']), ignore_index=True)
            
        elif element.element_type == "dipole" :
            B = element.Bfield
            L = element.length
            n = element.order
            k1 = element.k1
            k2 = element.k2
            gap = element.gap
            pole_face1 = element.angle1
            pole_face2 = element.angle2
            
            
            dipole_mat = sec_bend_matrix(L,B,n, self.get_p(), rest_mass=self.rest_mass)
            
            
            # add pole face calculations if non zero
            if ~np.isnan(pole_face1) :
                angle_mat1 = pole_face_mat(pole_face1, B, gap, self.get_p(), k1, k2)
                print(angle_mat1)
                dipole_mat = angle_mat1 * dipole_mat
            if ~np.isnan(pole_face2) :
                angle_mat2 = pole_face_mat(pole_face2, B, gap, self.get_p(), k1, k2)
                dipole_mat = dipole_mat * angle_mat2
                print(angle_mat2)
                
            # calculate new particle properties and put them in a new raw
            new_row = pd.DataFrame(data = [np.matmul(dipole_mat, self.get_vect())], \
                                     columns = self.p_df.columns)
            self.p_df = self.p_df.append(new_row, ignore_index=True)
            
            # update z
            self.z_df = self.z_df.append(pd.DataFrame(data = [self.get_z() + L], columns = ['z [m]']), ignore_index=True)
            
        
        elif element.element_type == "quad" :
            
            
        
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
        else:
            warnings.warn("Element transfer function not defined: ",element.element_type) 
            
            
            
    def particle_through_BL(self, BL : Beamline()):
        
        for index, row in BL.BL_df.iterrows(): # loop all elements into BL
            self.particle_through_BLelement(row['BL object'])
            


    


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
    


def create_BL_from_Transport(transport_file, scaling_ratio = 1):
    """
    opens a transport file and creates the corresponding Beamline object
    """
    
    # default settings
    my_beamline = Beamline()
    
    vertical_gap = 0.03 # default value
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       last_pos = fp.tell()
       
       prev_line = fp.readline()
       
       
       
       
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           data_prev = prev_line.split()
           
           
           if data: # string is not empty
               
               if data[0][0:2] == "3.":
                   # drift
                   L = float(data[1].replace(";","") )
                   
                   # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   new_element = BL_Element(name = name, length = L)
                   
                   my_beamline.add_element(new_element)
               
                    
               if data[0][0:2] == "4.":
                   # sector bending
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )*scaling_ratio
                   n = float(data[3].replace(";","") )
                   
                    # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   
                   angle_in = 0
                   angle_out = 0
                   if data_prev :
                       if data_prev[0] == "2.": # entrance pole face
                           angle_in = float(data_prev[1].replace(";","") )
                           
                   last_pos = fp.tell() # save line position
                   line = fp.readline()
                   if data_prev[0] == "2.": # exit pole face
                       angle_out = float(data_prev[1].replace(";","") )
                   else:
                       fp.seek(last_pos) # go back, exit face not defined
                       
                  
                   
                   new_element = BL_Element(name = name, length = L)
                   new_element.dipole(B=B, n=n, gap=vertical_gap, pole_face1=angle_in, pole_face2=angle_out)
                   my_beamline.add_element(new_element)
                
               if data[0][0:2] == "5.":
                   # quad
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") ) * scaling_ratio
                   a = float(data[3].replace(";","") ) / 1000 # aperture converted to meters
                   
                   # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   new_element = BL_Element(name = name, length = L)
                   new_element.quadrupole(B=B, a=a)
                   my_beamline.add_element(new_element)
                   
                   
                   
                   
               if data[0][0:3] == "16.":
                  if data[1][0:2] == "5." or data[1][0:1] == '5' :
                       # gap
                       vertical_gap = float(data[2].replace(";","") )/1000
                       
                       
           prev_line = line     
           line = fp.readline()
           
           
       return my_beamline

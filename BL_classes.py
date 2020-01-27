# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""


import numpy as np
import pandas as pd
from scipy.stats import norm
from transfer_functions import PtoGamma, PtoBrho, EtoP, PtoE, drift_matrix, quad_matrix, sec_bend_matrix, pole_face_mat, Highland
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt, pi
import warnings
import time
import copy

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
    
class BL_Element:
    """
    Define any element in the beamline.  Default = drift
    """

    def __init__(self, name="", length=0):
        self.name = name
        self.length = length
        self.element_type = "drift" # default is a drift
        self.N_segments = 1
        self.curvature = 0 # curvature of the element (cf. dipole case)
        self.aperture_type = "rectangular" # shape of element aperture (or vaccum chamber)
        self.apertureX = 1 # default aperture = 1m
        self.apertureY = 1 # default aperture = 1m
        
        
    def dipole(self, B, n, apertureY, apertureX=0, pole_face1=0, pole_face2=0, \
               curvature=0, CCT_angle = 0, k1 = 0.5, k2 = 0, \
               aperture_type = 'circular' ,nb_pts=10):
        self.element_type = "dipole"
        self.Bfield = B
        self.order = n
        self.apertureY = apertureY
        if apertureX == 0: self.apertureX = apertureY # unless specified, the aperture is symmetrical
        else: self.apertureX = apertureX
        self.aperture_type = aperture_type # rectanglar or circular
        self.k1 = k1
        self.k2 = k2
        self.pole_face1 = pole_face1
        self.pole_face2 = pole_face2
        self.N_segments = nb_pts
        self.curvature = curvature
        self.CCT_angle = CCT_angle # angle of the windings in the CCT case
        
    
    
    def set_curvature(B, p, rest_mass):
        """ determine curvature based on particle properties to the magnet field """
        return B * p * rest_mass
    
    def set_dipole_curvature(self, p, rest_mass):
        """ set dipole curvature based on particle properties to the magnet field """
        self.curvature = set_curvature(self.Bfield , p, rest_mass) 
        
    def quadrupole(self, B, a, CCT_angle = 0, nb_pts=10):
        self.element_type = "quad"
        self.Bfield = B
        self.apertureX = a
        self.apertureY = a
        self.aperture_type = "circular"
        self.N_segments = nb_pts
        self.CCT_angle = CCT_angle # angle of the windings in the CCT case
        
    def sextupole(self, B, a, CCT_angle = 0, nb_pts=10):
        self.element_type = "sext"
        self.Bfield = B
        self.apertureX = a
        self.apertureY = a
        self.aperture_type = "circular"
        self.N_segments = nb_pts
        self.CCT_angle = CCT_angle # angle of the windings in the CCT case
        
    def solenoid(self, B, a, nb_pts=20):
        self.element_type = "solenoid"
        self.Bfield = B
        self.apertureX = a
        self.apertureY = a
        self.aperture_type = "circular"
        self.N_segments = nb_pts
        
    def slit(self, gap, orientation='X', mat ='tantalum', offset = 0, nb_pts=20):
        self.element_type = "slit"
        
        if orientation == 'X':
            self.apertureX = gap/2
        elif orientation == 'Y':
            self.apertureY = gap/2
            
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
        slice_element.length = element.length
        
        
        new_row = pd.DataFrame(data = [[self.get_z() + slice_element.length, \
                                        slice_element.name, \
                                        slice_element.element_type, \
                                        slice_element ]], \
                                 columns = self.BL_df.columns)
        
        
        
        self.BL_df = self.BL_df.append(new_row, ignore_index=True)
    
    
    def get_element_index(self, element_name):
        """ returns the index of an element in the beamline """
        for index, row in self.BL_df[['Element name']].iterrows():   
            if row.values[0] == element_name:
                return index
            
        # out of loop --> element not found
        raise Exception("Element %s is not present in beamline"%element_name)
    
        

class Particle:
        
    def __init__(self, z=0, x=0, y=0, divX=0, divY=0, dL=0, p=0, refp=1, particle_type='proton', max_it=1000):
        """ 
        Create particle with some initial parameters 
        max_it = number of transformations the particle can go through in the beamline
        """
        
        if particle_type=='proton' : self.rest_mass = 938
        else: raise Exception("Particle type is not known")
        
        if refp <= 0:
            raise Exception("Reference energy is not a positive number")
        
        self.refp = refp
        self.particle_type = particle_type
        
        self.it = 0
        
        # initialize z coordinate
        self.z = np.empty(shape = [max_it, 1])
        self.z[0] = z
        
        # nitialize particle dynamic properties
        self.X = np.empty(shape = [max_it, 6])
        self.X[0, :] = [x, divX, y, divY, dL, (p-refp)/refp]
    
    def get_z(self, row_nb=-1): 
        if row_nb == -1:  return self.z[self.it] # last value that was updated
        else:   return self.z[row_nb]
    
    def get_x(self, row_nb=-1):
        if row_nb == -1:  return self.X[self.it, 0] # last value that was updated
        else:   return self.X[row_nb, 0]
        
    def get_y(self, row_nb=-1): 
        if row_nb == -1:  return self.X[self.it, 2] # last value that was updated
        else:   return self.X[row_nb, 2]
        
    def get_divx(self, row_nb=-1): 
        if row_nb == -1:  return self.X[self.it, 1] # last value that was updated
        else:   return self.X[row_nb, 1]
        
    def get_divy(self, row_nb=-1):
        if row_nb == -1:  return self.X[self.it, 3] # last value that was updated
        else:   return self.X[row_nb, 3]
        
    def get_dL(self, row_nb=-1): 
        if row_nb == -1:  return self.X[self.it, 4] # last value that was updated
        else:   return self.X[row_nb, 4]
        
    def get_dponp(self, row_nb=-1): 
        if row_nb == -1:  return self.X[self.it, 5] # last value that was updated
        else:   return self.X[row_nb, 5]
        
    def get_vect(self, row_nb=-1): 
        # get beam properties at row=row_nb, default=last row
        if row_nb == -1:  return self.X[self.it, :] # last value that was updated
        else:   return self.X[row_nb, :]
    
    def get_p(self, row_nb=-1): return self.get_dponp(row_nb)*self.refp + self.refp
    def get_E(self, row_nb=-1): return PtoE(self.get_p(row_nb))
    
    def plot_x(self): # example, to put in separate function later
        plt.plot(self.z, self.X[:, 0])
        
    
    def particle_through_drift(self, drift : BL_Element):
        L = drift.length
        drift_mat = drift_matrix(L, self.get_p(), rest_mass=self.rest_mass)
        
        # calculate new particle properties 
        self.it = self.it + 1
        self.X[self.it, :] = np.matmul(drift_mat, self.X[self.it - 1, :])
        self.z[self.it] = self.z[self.it-1] + L
        
        return 
    
    
    def particle_through_slit(self, slit : BL_Element):
        """
        slit with 2 blades (left and right)
        """
        
        L_segment = slit.length/slit.N_segments
        
        if slit.orientation == 'X':
            left_blade = -slit.apertureX/2 + slit.offset
            right_blade = slit.apertureX/2 + slit.offset
        elif slit.orientation == 'Y':
            left_blade = -slit.apertureY/2 + slit.offset
            right_blade = slit.apertureY/2 + slit.offset
        
        for i in range(0, slit.N_segments):
            if slit.orientation == 'X':
                pos = self.get_x()
            elif slit.orientation == 'Y':
                pos = self.get_y()
            
            
            if (pos <= left_blade or pos >= right_blade): # particle hits the slit
                    
                self.it = self.it + 1
                
                [self.X[self.it, :] , sigma] = Highland(L_segment, slit.material, self.X[self.it - 1, :], PtoE(self.refp))
                
                self.z[self.it] = self.z[self.it-1] + L_segment
                
                
            else: # particle passes through --> drift
                self.particle_through_drift(BL_Element(length=L_segment))
        
        return 
    
    
    
    def kill_lost_particle(self, element: BL_Element):
        if element.aperture_type == "rectangular":
            if abs(self.X[self.it, 0]) > element.apertureX or abs(self.X[self.it, 2]) > element.apertureY :
                # kill particle
                self.X[self.it, :] = np.nan
        elif element.aperture_type == "circular":
            if (self.X[self.it, 0])**2 + (self.X[self.it, 2])**2 > element.apertureY**2:
                # kill particle
                self.X[self.it, :] = np.nan
    
    def particle_through_BLelement(self, element : BL_Element):
        
        if np.isnan(self.get_dponp()) :
            # add rows with NAN
            self.it = self.it + 1
            self.X[self.it, :] = np.nan
            self.z[self.it] = self.z[self.it-1] + element.length
            return
        
        if element.element_type == "drift" :            
            self.particle_through_drift(element)
            
        elif element.element_type == "dipole" :
            
            N_segments = element.N_segments
            
            B = element.Bfield
            L = element.length/N_segments
            n = element.order
            k1 = element.k1
            k2 = element.k2
            apertureY = element.apertureY
            pole_face1 = element.pole_face1
            pole_face2 = element.pole_face2
            
            
            dipole_mat = sec_bend_matrix(L,B,n, self.get_p(), rest_mass=self.rest_mass)
            
            
            # entrance pole face calculations if non zero
            if ~np.isnan(pole_face1) :
                angle_mat1 = pole_face_mat(pole_face1, B, apertureY, self.get_p(), k1, k2)
                
                # calculate new particle properties 
                self.it = self.it + 1
                self.X[self.it, :] = np.matmul(angle_mat1, self.X[self.it - 1, :])
                self.z[self.it] = self.z[self.it-1] 
                
                
            
            # insisde dipole
            for i in np.arange(0, N_segments):
                # calculate new particle properties and put them in a new raw
                
                self.it = self.it + 1
                self.X[self.it, :] = np.matmul(dipole_mat, self.X[self.it - 1, :])
                self.z[self.it] = self.z[self.it-1] + L
                
                self.kill_lost_particle(element)
                
            # exit pole face calculations if non zero
            if ~np.isnan(pole_face2) :
                angle_mat2 = pole_face_mat(pole_face2, B, apertureY, self.get_p(), k1, k2)
                
                # calculate new particle properties 
                self.it = self.it + 1
                self.X[self.it, :] = np.matmul(angle_mat2, self.X[self.it - 1, :])
                self.z[self.it] = self.z[self.it-1] 
                
        elif element.element_type == "quad" :
            
            N_segments = element.N_segments
        
            B = element.Bfield
            L = element.length/N_segments
            a = element.apertureY
            
            quad_mat = quad_matrix(L, B, a, self.get_p(), rest_mass=self.rest_mass)
            
            
            for i in np.arange(0, N_segments):
                # calculate new particle properties 
                self.it = self.it + 1
                self.X[self.it, :] = np.matmul(quad_mat, self.X[self.it - 1, :])
                self.z[self.it] = self.z[self.it-1] + L
                
                self.kill_lost_particle(element)
                
        elif element.element_type == "slit" :                
            self.particle_through_slit(element)
        else:
            warnings.warn("Element transfer function not defined: ",element.element_type) 
            
        return    
            
    def particle_through_BL(self, BL : Beamline()):
        
        for index, row in BL.BL_df.iterrows(): # loop all elements into BL
            #print(" \n  BL element: %s (%s) ; z =  %.2e"%(row['Element name'],row['Element type'],BL.get_z(index)))
            
            self.particle_through_BLelement(row['BL object'])
            #print(self.get_z())


class Beam:
    
    def __init__(self, nb_part=1000, \
                       refE = 160, DeltaE=0, E_dist='cst',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform', \
                       particle_type='proton'):
        
        # initialize empty list of particles
        self.particle_list = []
        
        for i in range(0,nb_part):
            # initialize particle properties
            if size_dist =='uniform':
                sizeX = DeltaX*np.random.uniform(-1,1)
                sizeY = DeltaY*np.random.uniform(-1,1)
            elif size_dist == 'normal':
                sizeX = DeltaX*np.random.normal(0,1)
                sizeY = DeltaY*np.random.normal(0,1)
            elif size_dist == 'cst':
                sizeX = DeltaX*(i-nb_part/2)/nb_part*2
                sizeY = DeltaY*(i-nb_part/2)/nb_part*2
            else:
                raise Exception('Distribution chosen for "size_dist" is not valid')
            
            
            if div_dist =='uniform' :
                divX = DeltaDivX*np.random.uniform(-1,1)
                divY = DeltaDivY*np.random.uniform(-1,1)
            elif div_dist == 'normal':
                divX = DeltaDivX*np.random.normal(0,1)
                divY = DeltaDivY*np.random.normal(0,1)
            elif div_dist == 'cst':
                divX = DeltaDivX*(i-nb_part/2)/nb_part*2
                divY = DeltaDivY*(i-nb_part/2)/nb_part*2
            else:
                raise Exception('Distribution chosen for "div_dist" is not valid')
                
             
            if E_dist == 'uniform' :
                E = refE + DeltaE*np.random.uniform(-1,1)
            elif E_dist == 'normal':
                E = refE + DeltaE*np.random.normal(0,1)
            elif E_dist == 'cst':
                E = refE + DeltaE*(i-nb_part/2)/nb_part*2
            else:
                raise Exception('Distribution chosen for "E_dist" is not valid')
                
            
            # create new particle
            new_particle = Particle(z=0, x=sizeX , y=sizeY, divX=divX, divY=divY, \
                                    dL=0, p=EtoP(E), refp=EtoP(refE), particle_type=particle_type)
                                                                                                   
        
            self.particle_list = self.particle_list + [new_particle]
    
        
    def add_particle(self, particle : Particle()):
        self.particle_list = self.particle_list + [particle]
        
    def get_beam_X(self, row_nb=-1):
        """ get the X coordinates of all the particles in the beam in a given row"""
        
        X = [] # list of particle positions in X
        for particle in self.particle_list :
            X.append(particle.get_x(row_nb))
        
        return np.asarray(X) # convert to numpy
    
    
    def size_X(self, row_nb=-1):
        X = self.get_beam_X(row_nb)
        (muX, sigmaX) = norm.fit(X[~np.isnan(X)])
        return sigmaX
    
    


def create_BL_from_Transport(transport_file, scaling_ratio = 1, CCT_angle=0):
    """
    opens a transport file and creates the corresponding Beamline object
    """
    
    # default settings
    my_beamline = Beamline()
    
    vertical_aperture = 0.03 # default value
    
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
               
                    
               elif data[0][0:2] == "4.":
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
                   data = line.split() 
                   
                   if data: # string is not empty
                       if data_prev[0] == "2.": # exit pole face
                           angle_out = float(data_prev[1].replace(";","") )
                       else:
                           fp.seek(last_pos) # go back, exit face not defined
                       
                  
                   
                   new_element = BL_Element(name = name, length = L)
                   new_element.dipole(B=B, n=n, apertureY=vertical_aperture, \
                                      pole_face1=angle_in, pole_face2=angle_out,\
                                      CCT_angle = CCT_angle)
                   my_beamline.add_element(new_element)
                
               elif data[0][0:2] == "5.":
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
                   new_element.quadrupole(B=B, a=a, CCT_angle = CCT_angle)
                   my_beamline.add_element(new_element)
                   
               elif data[0][0:5] == '(slit':
                   # slit
                   
                   L = float(data[1].replace(";","") )
                   leftX = float(data[2].replace(";","") )
                   rightX = float(data[3].replace(";","") )
                   mat = data[4].replace(";","") 
                   
                   opening = rightX - leftX
                   offset = rightX + leftX
                   
                   new_element = BL_Element(name = name, length = L)
                   new_element.slit(opening, orientation='X', mat = mat, offset = offset, nb_pts=20)
                   my_beamline.add_element(new_element)
                   
                   line = fp.readline() # skip next drift
                   
                   
               elif data[0][0:3] == "16.":
                  if data[1][0:2] == "5." or data[1][0:1] == '5' :
                       # aperture
                       vertical_aperture = float(data[2].replace(";","") )/1000
                   
                    
                    
                       
           prev_line = line     
           line = fp.readline()
           
           
       return my_beamline




    






def magnet_patches(element: BL_Element, orientation = 'ZX', coil_height=0, start_point = np.array([0,0])):
    """ 
    Returns coordinates of element used to plot the magnet
    Possible orientations for plot: ZX or ZY
    """
    
    x0 = start_point[0]
    y0 = start_point[1]
    
    
    if orientation == 'ZX':
        if hasattr(element, 'apertureX'):
            a = element.apertureX
        else:
            return []
    elif orientation == 'ZY':
        if hasattr(element, 'apertureY'):
            a = element.apertureY
        else:
            return []
    
    
    
    if coil_height > 0:
        h = coil_height + 0.02 # add 2 cm for formers/insulation/etc
    else:
        h = 2 * a + 0.02 # default value for coil representation ; add 2 cm for formers
    L = element.length
    
    if hasattr(element, 'CCT_angle'):
        if element.element_type == 'dipole':
            CCT_angle = element.CCT_angle
        elif  element.element_type == 'quad':
            CCT_angle = element.CCT_angle/2
        elif  element.element_type == 'sext':
            CCT_angle = element.CCT_angle/3    
    else:
        CCT_angle = 0
    
    upper_points = [[x0, y0 + a] , [x0 + L, y0 + a] , [x0 + L + tan(CCT_angle)*h, y0 + a + h], [x0 - tan(CCT_angle)*h, y0 + a + h] ]
    lower_points = [[x0, y0 - a] , [x0 + L, y0 - a] , [x0 + L + tan(CCT_angle)*h, y0 - a - h], [x0 - tan(CCT_angle)*h, y0 - a - h] ]
    
    patches = []
    
    polygon1 = Polygon(np.array(upper_points), closed = True)
    patches.append(polygon1)
    
    polygon2 = Polygon(np.array(lower_points), closed = True)
    patches.append(polygon2)
    
    # coordinates for labelling
    center_x = x0 + L/2
    center_y = y0 + a + h/2
    
    return [patches, center_x, center_y]
        
def BL_plot(BL : Beamline):
    """ Possible orientations for plot: ZX or ZY """
    
    start_point = np.array([0,0])
    
    element_list = BL.BL_df['BL object'].values
    
    fig = plt.figure('Beamline',figsize=(18, 12))
    
    # increase vertical space between plots
    fig.subplots_adjust(hspace=0.4)
    
    ax_X = fig.add_subplot(2,1,1)
    ax_Y = fig.add_subplot(2,1,2)
    
    for element in element_list:
        
        if element.element_type == 'dipole' :
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='blue', alpha=0.6)
            ax_X.add_collection(p)
            
            ax_X.text(center_x, center_y, element.name, \
                      verticalalignment='center', horizontalalignment='center')
            
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='blue', alpha=0.6)
            ax_Y.add_collection(p)
            ax_Y.text(center_x, center_y, element.name, \
                      verticalalignment='center', horizontalalignment='center')
            
        elif element.element_type == 'quad' :
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='red', alpha=0.6)
            ax_X.add_collection(p) 
            ax_X.text(center_x, center_y, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
            
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='red', alpha=0.6)
            ax_Y.add_collection(p)
            ax_Y.text(center_x, center_y, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
        elif element.element_type == 'slit' :
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='dimgray', alpha=0.8)
            ax_X.add_collection(p)   
            
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='dimgray', alpha=0.8)
            ax_Y.add_collection(p)   
        
            
        start_point = start_point + [element.length , 0]
        
    ax_X.set_xlim([-0.1, start_point[0]])
    ax_X.set_ylim([-0.15,0.15])
    ax_X.set_xlabel('z [m]')
    ax_X.set_ylabel('x [m]')
    
    ax_Y.set_xlim([-0.1, start_point[0]])
    ax_Y.set_ylim([-0.15,0.15])
    ax_Y.set_xlabel('z [m]')
    ax_Y.set_ylabel('y [m]')
    
    return [fig, ax_X, ax_Y]


        
        
        
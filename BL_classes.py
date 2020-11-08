# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:40:04 2020

@author: rvroerm
"""


import numpy as np
import pandas as pd
from scipy.stats import norm
from transfer_functions import PtoGamma, PtoBrho, EtoP, PtoE, drift_matrix, quad_matrix, sec_bend_matrix, pole_face_mat, solenoid_mat, Highland, rot_mat
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt, pi, atan, asin, acos
import warnings
import time
import copy

import matplotlib.pyplot as plt


    
class BL_Element:
    """
    Define beamline element.  Default = drift, other elements in child classes
    """

    def __init__(self, name="", length=0, element_rot_rad=0):
        self.name = name
        self.length = length
        self.element_type = "drift" # default is a drift
        self.N_segments = 1
        self.curvature = 0 # curvature of the element (cf. dipole case)
        self.aperture_type = "rectangular" # shape of element aperture (or vaccum chamber)
        self.apertureX = 1 # default aperture = 1m
        self.apertureY = 1 # default aperture = 1m
        self.coil_height = 0 # default coil height
        self.element_rot_rad = element_rot_rad # rotation of the element in the X/Y plane
     
    
    def set_coil_height(self, h):
        self.coil_height = h
    
    def set_curvature(self, B, p):
        """ determine curvature based on particle properties to the magnet field """
        self.curvature = B / (p/300) # assumes p in Mev/c
        
    def get_bending_angle(self):
        if self.curvature == 0: warnings.warn("Magnet curvature doesn't seem to be set")
        return self.curvature * self.length
    
        
    
class Dipole(BL_Element):
    """
    Dipole class
    """
    
    def __init__(self, name="", length=0, element_rot_rad=0, B=0, n=0, apertureY=0, apertureX=0, pole_face1=0, pole_face2=0, \
               curvature=0, CCT_angle = 0, k1 = 0.5, k2 = 0, \
               aperture_type = 'circular' ,nb_pts=10):
        super().__init__(name=name, length=length, element_rot_rad=element_rot_rad)
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
		
		
class Quadrupole(BL_Element):
    
    def __init__(self, name="", length=0, element_rot_rad=0, B=0, a=0, CCT_angle = 0, \
               aperture_type = 'circular' ,nb_pts=10):
        super().__init__(name=name, length=length, element_rot_rad=element_rot_rad)
        self.element_type = "quad"
        self.Bfield = B
        self.apertureY = a
        self.apertureX = a
        self.aperture_type = aperture_type
        self.N_segments = nb_pts
        self.CCT_angle = CCT_angle # angle of the windings in the CCT case
        
class Sextupole(BL_Element):
    
    def __init__(self, name="", length=0, element_rot_rad=0, B=0, a=0, CCT_angle = 0, \
               aperture_type = 'circular' ,nb_pts=10):
        super().__init__(name=name, length=length, element_rot_rad=element_rot_rad)
        self.element_type = "sext"
        self.Bfield = B
        self.apertureY = a
        self.apertureX = a
        self.aperture_type = aperture_type
        self.N_segments = nb_pts
        self.CCT_angle = CCT_angle # angle of the windings in the CCT case

class ScanMagnet(BL_Element):
    
    def __init__(self, name="", length=0, Bx=0, By=0, apertureX=0, apertureY=0, \
               aperture_type = 'rectangular' ,nb_pts=10):
        super().__init__(name=name, length=length)
        self.element_type = "SM"
        self.Bx = Bx
        self.By = By
        self.apertureX = apertureX
        self.apertureY = apertureY
        self.aperture_type = aperture_type
        self.N_segments = nb_pts
        
class Solenoid(BL_Element):
    
    def __init__(self, name="", length=0, B=0, a=0, \
               aperture_type = 'circular' ,nb_pts=20):
        super().__init__(name=name, length=length)
        self.element_type = "solenoid"
        self.Bfield = B
        self.apertureY = a
        self.apertureX = a
        self.aperture_type = aperture_type
        self.N_segments = nb_pts
        
        
class Slit(BL_Element):
    
    def __init__(self, gap=0, name="", length=0, element_rot_rad=0, orientation='X', mat ='lead', offset = 0, nb_pts=20):
        
        super().__init__(name=name, length=length, element_rot_rad=element_rot_rad)
        self.element_type = "slit"
        
        if orientation == 'X':
            self.apertureX = gap/2
        elif orientation == 'Y':
            self.apertureY = gap/2
        self.orientation = orientation
        self.material = mat
        self.offset = offset
        self.N_segments = nb_pts
        
class BPM(BL_Element):
    
    def __init__(self, name="", length=0):
        super().__init__(name=name, length=length)
        self.element_type = "BPM"
        
        
##############################################################################        

class Beamline:
    
    def __init__(self):
        
        
        self.BL_df = pd.DataFrame([[ 0, "souce", "", BL_Element()]], \
                                  columns = ['z [m]', 'Element name', 'Element type', 'BL object'])
        
        
        
        
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
    
    def get_total_segments(self):
        N = 1
        for index, row in self.BL_df.iterrows(): # loop all elements into BL
            
            N = N + row['BL object'].N_segments
            if row['BL object'].element_type == 'dipole':
                N = N + 2 # add 2 segements for the pole faces
        
        return N
    
    



###############################################################################


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
        self.z = np.empty(shape = [max_it,])
        self.z[:] = np.nan
        self.z[0] = z
        
        # nitialize particle dynamic properties
        self.X = np.empty(shape = [max_it, 6])
        self.X[:, :] = np.nan
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
    
    def get_param(self, param='x', row_nb=-1):
        if param == 'z' : return self.get_z(row_nb)
        elif param == 'x' : return self.get_x(row_nb)
        elif param == 'y' : return self.get_y(row_nb)
        elif param == 'divx' : return self.get_divx(row_nb)
        elif param == 'divy' : return self.get_divy(row_nb)
        elif param == 'dL' : return self.get_dL(row_nb)
        elif param == 'dponp' : return self.get_dponp(row_nb)
        elif param == 'p' : return self.get_p(row_nb)
        elif param == 'E' : return self.get_E(row_nb)
        else: raise Exception("Parameter %s is not valid"%param) 
    
    def get_z_index(self, z_lookup, tolerance =10**-6):
        """ return index of a given value of z (z_lookup) """
        
        # remove nan from z list
        z_list = self.z[:] 
        z_list[np.isnan(z_list)] = 0
        
        p_index = np.where( abs(z_list - z_lookup) < tolerance )
        if len(p_index[0]) == 1 :
            p_index = np.asscalar(p_index[0])
        else:
            raise Exception('ISO not found')
            
        return p_index
    
    def particle_through_drift(self, drift : BL_Element):
        L = drift.length
        drift_mat = drift_matrix(L, self.get_p(), rest_mass=self.rest_mass)
        
        # calculate new particle properties 
        self.it = self.it + 1
        self.X[self.it, :] = np.matmul(drift_mat, self.X[self.it - 1, :])
        self.z[self.it] = self.z[self.it-1] + L
        
        return 
    
    
    def particle_through_slit(self, slit : Slit):
        """
        slit with 2 blades (left and right)
        """
        
        L_segment = slit.length/slit.N_segments
        
        if slit.orientation == 'X':
            left_blade = -slit.apertureX - slit.offset
            right_blade = slit.apertureX - slit.offset
        elif slit.orientation == 'Y':
            left_blade = -slit.apertureY - slit.offset
            right_blade = slit.apertureY - slit.offset
        
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
    
    
        
    
    def particle_through_dipole(self, element: Dipole):
        N_segments = element.N_segments 
            
        B = element.Bfield
        L = element.length/N_segments
        n = element.order
        k1 = element.k1
        k2 = element.k2
        apertureY = element.apertureY
        pole_face1 = element.pole_face1
        pole_face2 = element.pole_face2
        
        # X/Y rotation of element
        angle = element.element_rot_rad
        rot_mat1 = rot_mat(angle)
        rot_mat2 = rot_mat(-angle)
        
        dipole_mat = rot_mat2 @ sec_bend_matrix(L,B,n, self.get_p(), rest_mass=self.rest_mass) @ rot_mat1
        
        # entrance pole face calculations if non zero
        if ~np.isnan(pole_face1) :
            angle_mat1 = rot_mat2 @ pole_face_mat(pole_face1, B, apertureY, self.get_p(), k1, k2) @ rot_mat1
            
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
            angle_mat2 = rot_mat2 @ pole_face_mat(pole_face2, B, apertureY, self.get_p(), k1, k2) @ rot_mat1
            
            # calculate new particle properties 
            self.it = self.it + 1
            self.X[self.it, :] = np.matmul(angle_mat2, self.X[self.it - 1, :])
            self.z[self.it] = self.z[self.it-1] 
            
        return
    
    
    def particle_through_quadrupole(self, element: Quadrupole):
        N_segments = element.N_segments
        
        B = element.Bfield
        L = element.length/N_segments
        a = element.apertureY
        
        # X/Y rotation of element
        angle = element.element_rot_rad
        rot_mat1 = rot_mat(angle)
        rot_mat2 = rot_mat(-angle)
        
        quad_mat = rot_mat2 @ quad_matrix(L, B, a, self.get_p(), rest_mass=self.rest_mass) @ rot_mat1
        
        
        for i in np.arange(0, N_segments):
            # calculate new particle properties 
            self.it = self.it + 1
            self.X[self.it, :] = np.matmul(quad_mat, self.X[self.it - 1, :])
            self.z[self.it] = self.z[self.it-1] + L
            
            self.kill_lost_particle(element)
            
        return
    
    def particle_through_solenoid(self, element: Solenoid):
        N_segments = element.N_segments
        
        B = element.Bfield
        L = element.length/N_segments
        a = element.apertureX
        
        solen_mat = solenoid_mat(L, B, a,beam,it_z,refE, N_segments=N_segments)
        
        for i in np.arange(0, N_segments):
            # calculate new particle properties 
            self.it = self.it + 1
            self.X[self.it, :] = np.matmul(solenoid_mat, self.X[self.it - 1, :])
            self.z[self.it] = self.z[self.it-1] + L
            
            self.kill_lost_particle(element)
    
    def kill_lost_particle(self, element: BL_Element):
        if element.aperture_type == "rectangular":
            if abs(self.get_x(self.it)) > element.apertureX or abs(self.get_y(self.it)) > element.apertureY :
                # kill particle
                self.X[self.it, :] = np.nan
        elif element.aperture_type == "circular":
            if (self.get_x(self.it))**2 + (self.get_y(self.it))**2 > element.apertureY**2:
                # kill particle
                self.X[self.it, :] = np.nan
    
    
    def particle_through_SM(self, element: ScanMagnet):
        " test scanning magnet - rectangular"
        
        L_segment = element.length / element.N_segments
        p = self.get_p()
        
        Bx = element.Bx
        By = element.By
        Bz = 0 # assume transverse field
        B = sqrt(Bx**2 + By**2 + Bz**2)
        
        if B==0:
            self.particle_through_drift(BL_Element(length=element.length))
            return
        
        # get parallel and orthogonal components of p vs B
        hypotenuse_p = sqrt(self.get_divx()**2 + self.get_divy()**2 + 1)
        px = self.get_divx()/hypotenuse_p * p
        py = self.get_divy()/hypotenuse_p * p
        pz = 1/hypotenuse_p * p
       
        p_parallel = (px*Bx + py*By + pz*Bz) / B
        p_orth = sqrt(p**2 - p_parallel**2)
        
        # get curvature radius
        Brho  = PtoBrho(p_orth)
        rho = Brho/B 
        
        
        # rotate B on +y axis
        if By == 0:
            By_angle = pi/2 + pi*(Bx<0) 
        else:
            By_angle = atan(Bx/By) - pi*(By>0)
        
        rot_Byangle = np.array([[cos(By_angle) , -sin(By_angle)],
                                [sin(By_angle) , cos(By_angle)]])
        
        rot_Byangle2 = np.array([[cos(By_angle) , sin(By_angle)],
                                [-sin(By_angle) , cos(By_angle)]])
        
        # rotate particle initial position and angle in rotated frame
        new_pos_vect = np.matmul(rot_Byangle, 
                                 np.array([[ self.get_x() ], [ self.get_y() ]]))  
        x = new_pos_vect[0]
        y = new_pos_vect[1]
        new_angle_vect = np.matmul(rot_Byangle, 
                                   np.array([[ self.get_divx() ], [ self.get_divy() ]]))  
        divX = new_angle_vect[0]
        divY = new_angle_vect[1]
        
        
        
        for i in range(0, element.N_segments):
            # get exit position and angle
            x = x + rho*cos(divX) - sqrt( rho**2 - (L_segment + rho*sin(divX))**2 )
            y = y + divY * L_segment # drift in direction parallel to B
            theta = 2 * asin(1/2/rho * sqrt( x**2 + L_segment**2 ))
            divX = divX + theta
            
            
            # rotate back
            pos_vect = np.matmul(rot_Byangle2, np.array([x, y]))
            div_vect = np.matmul(rot_Byangle2, np.array([divX, divY]))
            
            # update matrix
            self.it = self.it + 1
            self.X[self.it, 0] = pos_vect[0]
            self.X[self.it, 1] = div_vect[0]
            self.X[self.it, 2] = pos_vect[1]
            self.X[self.it, 3] = div_vect[1]
            self.X[self.it, 4] = self.X[self.it-1, 4]
            self.X[self.it, 5] = self.X[self.it-1, 5]
            self.z[self.it] = self.z[self.it-1] + L_segment
            
            
            self.kill_lost_particle(element)
        
        return
        
        
    def particle_through_BLelement(self, element : BL_Element):
        """
        computes the tranport of the particle through an beamline element 
        (drift, dipole, quadrupole, etc.)
        """
        
        if element.element_type == "drift" or element.element_type == "BPM" :            
            self.particle_through_drift(element)
            
        elif element.element_type == "dipole" :
            self.particle_through_dipole(element)
             
        elif element.element_type == "quad" :
            self.particle_through_quadrupole(element)
        
        elif element.element_type == "SM" :  
            self.particle_through_SM(element)
            
        elif element.element_type == "slit" :                
            self.particle_through_slit(element)
            
        else:
            warnings.warn("Element transfer function not defined: ",element.element_type) 
            
        return    
        
    
    def particle_through_BL(self, BL : Beamline()):
        """
        computes the tranport of the particle through all element in the beamline
        """
        
        if len(self.z) < BL.get_total_segments() :
            raise Exception("Beamline requires bigger length for particle set: increase max_it")
            
        for index, row in BL.BL_df.iterrows(): # loop all elements into BL
            #print(" \n  BL element: %s (%s) ; z =  %.2e"%(row['Element name'],row['Element type'],BL.get_z(index)))
            
            self.particle_through_BLelement(row['BL object'])
            #print(self.get_z())



###############################################################################

class Beam:
    
    
    
    def __init__(self, nb_part=1000, \
                       refE = 160, DeltaE=0, E_dist='uniform',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='uniform', \
                       OffsetX = 0, OffsetY=0, \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform', \
                       particle_type='proton'):
        """
        Initializes a beam with a set of nb_part particles
        
        The parameters for these particles are set statistically according to the inputs
        
        refE = central energy
        DeltaE = energy span around central energy
        E_dist = energy distribution with parameters refE, DeltaE. Possible distributions are cst, uniform, uniform2 (=uniform without random), normal
        DeltaX/DeltaY = spread size in X/Y
        OffsetX/Y = position offset in X/Y
        size_dist = position distribution with parameters refE, DeltaE. Possible distributions are cst, uniform, uniform2 (=uniform without random), normal
        DeltaDivX/Y = spread in divergence in X/Y
        div_dist = divergence distribution with parameters refE, DeltaE. Possible distributions are cst, uniform, uniform2 (=uniform without random), normal
        particle_type
        """
        
        # initialize empty list of particles
        self.particle_list = []
        
        for i in range(0,nb_part):
            # initialize particle properties
            if size_dist =='cst':
                posX = DeltaX + OffsetX
                posY = DeltaY + OffsetY
            elif size_dist == 'normal':
                posX = DeltaX*np.random.normal(0,1) + OffsetX
                posY = DeltaY*np.random.normal(0,1) + OffsetY
            elif size_dist == 'uniform':
                posX = DeltaX*np.random.uniform(-1,1) + OffsetX
                posY = DeltaY*np.random.uniform(-1,1) + OffsetY
            elif size_dist == 'uniform2': 
                #uniform distribution without randomness (to use on 1 variable max)
                posX = DeltaX*(i-nb_part/2)/nb_part*2 + OffsetX
                posY = DeltaY*(i-nb_part/2)/nb_part*2 + OffsetY
            else:
                raise Exception('Distribution chosen for "size_dist" is not valid')
            
            
            if div_dist =='cst' :
                divX = DeltaDivX
                divY = DeltaDivY
            elif div_dist == 'normal':
                divX = DeltaDivX*np.random.normal(0,1)
                divY = DeltaDivY*np.random.normal(0,1)
            elif div_dist == 'uniform':
                divX = DeltaDivX*np.random.uniform(-1,1)
                divY = DeltaDivY*np.random.uniform(-1,1)
            elif div_dist == 'uniform2': 
                #uniform distribution without randomness (to use on 1 variable max)
                divX = DeltaDivX*(i-nb_part/2)/nb_part*2
                divY = DeltaDivY*(i-nb_part/2)/nb_part*2
            else:
                raise Exception('Distribution chosen for "div_dist" is not valid')
                
             
            if E_dist == 'cst' :
                E = refE + DeltaE
            elif E_dist == 'normal':
                E = refE + DeltaE*np.random.normal(0,1)
            elif E_dist == 'uniform':
                E = refE + DeltaE*np.random.uniform(-1,1)
            elif E_dist == 'uniform2':
                E = refE + DeltaE*(i-nb_part/2)/nb_part*2
            else:
                raise Exception('Distribution chosen for "E_dist" is not valid')
                
            
            # create new particle
            new_particle = Particle(z=0, x=posX , y=posY, divX=divX, divY=divY, \
                                    dL=0, p=EtoP(E), refp=EtoP(refE), particle_type=particle_type)
                                                                                                   
        
            self.particle_list = self.particle_list + [new_particle]
    
        
    def add_particle(self, particle : Particle()):
        self.particle_list = self.particle_list + [particle]
        
    def get_beam_x(self, row_nb=-1):
        """ get the X coordinates of all the particles in the beam in a given row"""
        
        x = [] # list of particle positions in X
        for particle in self.particle_list :
            x.append(particle.get_x(row_nb))
        
        return np.asarray(x) # convert to numpy
    
    def get_beam_param(self, param='x', row_nb=-1):
        """ 
        get the values for the parameter 'param' for all the particles in the beam in a given row
        possible values for param: (z, x, y, divx, divy, dL, dponp, p, E)
        """
        my_list = [] # list of values for chosen parameter
        for particle in self.particle_list :
            my_list.append(particle.get_param(param=param, row_nb=row_nb))
        
        return np.asarray(my_list) # convert to numpy
        
        
    def pos_X(self, row_nb=-1, clip_dist=1):
        """
    
        Parameters
        ----------
        row_nb : int, optional
            Row number (in Particle dataframe) at which the property is calculated. The default is -1.
        clip_dist : float, optional
            Use  particles within [-clip_dist,clip_dist] only. The default is 1.

        Returns
        -------
        mu : float
            Transverse position X of the beam.
            
        """
        
        values = self.get_beam_param(param='x', row_nb=row_nb)
        (mu, sigma) = norm.fit(np.clip(values[~np.isnan(values)], -clip_dist, clip_dist))
        
        
        return mu
    
    def size_X(self, row_nb=-1, clip_dist=1):
        """
    
        Parameters
        ----------
        row_nb : int, optional
            Row number (in Particle dataframe) at which the property is calculated. The default is -1.
        clip_dist : float, optional
            Use  particles within [-clip_dist,clip_dist] only. The default is 1.

        Returns
        -------
        sigma : float
            Transverse size X of the beam.
            
        """
        
        values = self.get_beam_param(param='x', row_nb=row_nb)
        (mu, sigma) = norm.fit(np.clip(values[~np.isnan(values)], -clip_dist, clip_dist))
        
        
        return sigma
    
    
    def pos_Y(self, row_nb=-1, clip_dist=1):
        """
    
        Parameters
        ----------
        row_nb : int, optional
            Row number (in Particle dataframe) at which the property is calculated. The default is -1.
        clip_dist : float, optional
            Use  particles within [-clip_dist,clip_dist] only. The default is 1.

        Returns
        -------
        sigma : float
            Transverse position Y of the beam.
            
        """
        
        values = self.get_beam_param(param='y', row_nb=row_nb)
        (mu, sigma) = norm.fit(np.clip(values[~np.isnan(values)], -clip_dist, clip_dist))
        return mu
    
    
    def size_Y(self, row_nb=-1, clip_dist=1):
        """
    
        Parameters
        ----------
        row_nb : int, optional
            Row number (in Particle dataframe) at which the property is calculated. The default is -1.
        clip_dist : float, optional
            Use  particles within [-clip_dist,clip_dist] only. The default is 1.

        Returns
        -------
        sigma : float
            Transverse size Y of the beam.
            
        """
        
        values = self.get_beam_param(param='y', row_nb=row_nb)
        (mu, sigma) = norm.fit(np.clip(values[~np.isnan(values)], -clip_dist, clip_dist))
        return sigma
    
    def beam_through_BL(self, BL : Beamline()):
        for particle in self.particle_list :
            particle.particle_through_BL(BL)




###############################################################################

def create_BL_from_Transport(transport_file, scaling_ratio = 1, CCT_angle=0, apertureX=0):
    """
    opens a Transport input file and creates the corresponding Beamline object
    
    Optional parameters:
        scaling ratio: multiplies all magnetic fields
        CCT_angle: considers a CCT magnet with a wire inclination equal to the angle
        aperture_X: can define the vertical aperture of the dipoles as an input
    """
    
    # default settings
    my_beamline = Beamline()
    
    vertical_aperture = 0.03 # default value
    element_angle = 0
    
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
                   if posL != -1 and posR != -1 : # name found
                       name = line[posL+1 : posR].strip()
                       if name[0:3] == "BPM" or name[0:2] == "IC" :
                           # dosimetry element (BPM)
                           new_element = BPM(name = name, length = L)
                       else: new_element = BL_Element(name = name, length = L)
                   else: new_element = BL_Element(name = "", length = L)
                   
                   my_beamline.add_element(new_element)
               
                    
               elif data[0][0:2] == "4.":
                   # sector bending
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )*scaling_ratio
                   n = float(data[3].replace(";","") )
                   
                    # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   if posL != -1 and posR != -1 : # name found
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
                       if data[0] == "2.": # exit pole face
                           angle_out = float(data[1].replace(";","") )
                       else:
                           fp.seek(last_pos) # go back, exit face not defined
                   
                   # set horizontal aperture of the dipole
                   if apertureX!=0: #value was given in input
                       horizontal_aperture = apertureX 
                       aperture_type = "rectangular"
                   elif CCT_angle!=0 : # CCT magnet --> symmetric aperture
                       horizontal_aperture = 0 
                       aperture_type = "circular"
                   else: # default value = 5cm (x2)
                       horizontal_aperture = 0.05 
                       aperture_type = "rectangular"
                   
                   
                   new_element = Dipole(name = name, length = L, B=B, n=n, \
                                      apertureX=horizontal_aperture, apertureY=vertical_aperture, \
                                      aperture_type = aperture_type, \
                                      pole_face1=angle_in, pole_face2=angle_out,\
                                      CCT_angle = CCT_angle, element_rot_rad=element_angle)
                       
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
                   
                   new_element = Quadrupole(name = name, length = L, B=B, a=a, CCT_angle = CCT_angle, element_rot_rad=element_angle)
                   
                   my_beamline.add_element(new_element)
               
                
               elif data[0][0:3] == "19.":
                   # solenoid
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") ) * scaling_ratio
                   
                   if apertureX!=0: #value was given in input
                       a = apertureX 
                   elif vertical_aperture != 0:
                       a = vertical_aperture
                   else:
                       a = apertureY
                   
                   # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   new_element = Solenoid(name = name, length = L, B=B, a=a)
                   my_beamline.add_element(new_element)
                   
                
               elif data[0][0:3] == '(SM':
                   # scanning magnet
                   
                   L = float(data[1].replace(";","") )
                   Bx = float(data[2].replace(";","") ) * scaling_ratio
                   By = float(data[3].replace(";","") ) * scaling_ratio
                   apertureX = float(data[4].replace(";","") )
                   apertureY = float(data[5].replace(";","") )
                   
                   # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   new_element = ScanMagnet(name = name, length = L, Bx=Bx, By=By, apertureX=apertureX, apertureY=apertureY)
                   my_beamline.add_element(new_element)
                   
                   line = fp.readline() # skip next drift
                   
               elif data[0][0:5] == '(slit' or data[0][0:6] == '(slitX' or data[0][0:6] == '(slitY':
                   # slit
                   oritentation = 'X'
                   if data[0][0:6] == '(slitY':
                       oritentation = 'Y'
                   
                   L = float(data[1].replace(";","") )
                   left = float(data[2].replace(";","") )
                   right = float(data[3].replace(";","") )
                   mat = data[4].replace(";","") 
                   
                   opening = right - left
                   offset = right + left
                   
                   # check if there is a name
                   posL = line.find('/')
                   posR = line.rfind('/')
                   
                   if posL != -1 and posR != -1 :
                       name = line[posL+1 : posR].strip()
                   else: name = ""
                   
                   new_element = Slit(name = name, length = L, gap=opening, orientation=oritentation, mat = mat, offset = offset, nb_pts=20, element_rot_rad=element_angle)
                   my_beamline.add_element(new_element)
                   
                   line = fp.readline() # skip next drift
                   
                   
               elif data[0][0:3] == "16.":
                  if data[1][0:2] == "5." or data[1][0:1] == '5' :
                       # aperture
                       vertical_aperture = float(data[2].replace(";","") )/1000
                       
                       
               elif data[0][0:3] == "20.":
                   element_angle = element_angle + np.radians(float(data[1].replace(";","") ))
                   
                   
                   
                    
                    
                       
           prev_line = line     
           line = fp.readline()
           
           
       return my_beamline




    







        



        
        
        
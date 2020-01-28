# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:39:31 2020

@author: rvroerm
"""

from BL_classes import Particle, Beam, Beamline, BL_Element
from math import pi, tan
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import seaborn as sns
import time
import numpy as np
from scipy.stats import norm
from transfer_functions import EtoP

plt.close('all')


def plot_beam_through_BL(my_beamline: Beamline, my_beam: Beam):
    """ plot beam through BL """
    
    [fig, ax_X, ax_Y] = BL_plot(my_beamline)
    
    for particle in my_beam.particle_list :
        
        particle.particle_through_BL(my_beamline)
        
        ax_X.plot(particle.z[0:particle.it], particle.X[0:particle.it,0])
        ax_Y.plot(particle.z[0:particle.it], particle.X[0:particle.it,2])
        
    return [fig, ax_X, ax_Y]


def BL_plot(BL : Beamline):
    """ Possible orientations for plot: ZX or ZY """
    
    start_point = np.array([0,0])
    
    element_list = BL.BL_df['BL object'].values
    
    fig = plt.figure('Beamline',figsize=(18, 12))
    
    # increase vertical space between plots
    fig.subplots_adjust(hspace=0.4)
    
    ax_X = fig.add_subplot(2,1,1)
    ax_Y = fig.add_subplot(2,1,2)
    ylim = 0.15
    
    for element in element_list:
        
        if element.element_type == 'drift' and element.name!="" :
            ax_X.text(start_point[0] + element.length, ylim*0.8, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
            ax_X.axvline(x=start_point[0] + element.length, ymin=0.25, ymax=0.75, linestyle='-.', color='grey')
            
            ax_Y.text(start_point[0] + element.length, ylim*0.8, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
            ax_Y.axvline(x=start_point[0] + element.length, ymin=0.25, ymax=0.75, linestyle='-.', color='grey')
            
        elif element.element_type == 'dipole' :
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
            ax_X.text(center_x, ylim*0.8, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
            
            [patches, center_x, center_y] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='dimgray', alpha=0.8)
            ax_Y.add_collection(p) 
            ax_Y.text(center_x, ylim*0.8, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
        
            
        start_point = start_point + [element.length , 0]
        
    ax_X.set_xlim([-0.1, start_point[0]])
    ax_X.set_ylim([-ylim,ylim])
    ax_X.set_xlabel('z [m]')
    ax_X.set_ylabel('x [m]')
    
    ax_Y.set_xlim([-0.1, start_point[0]])
    ax_Y.set_ylim([-ylim,ylim])
    ax_Y.set_xlabel('z [m]')
    ax_Y.set_ylabel('y [m]')
    
    return [fig, ax_X, ax_Y]


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
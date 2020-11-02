# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:39:31 2020

@author: rvroerm
"""

from BL_classes import Beam, Beamline, BL_Element
from math import pi, sin, cos, tan
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection
#from planar import Polygon


import seaborn as sns
import numpy as np
from scipy.stats import norm
from transfer_functions import EtoP
import warnings

plt.close('all')


def plot_beam_through_BL(my_beamline: Beamline, my_beam: Beam):
    """ plot beam through BL """
    
    [fig, ax_X, ax_Y] = BL_plot_for_traces(my_beamline)
    
    for particle in my_beam.particle_list :
        
        particle.particle_through_BL(my_beamline)
        
        ax_X.plot(particle.z[0:particle.it], particle.X[0:particle.it,0])
        ax_Y.plot(particle.z[0:particle.it], particle.X[0:particle.it,2])
        
    return [fig, ax_X, ax_Y]


def BL_plot_for_traces(BL : Beamline, title='Beamline'):
    """ 
    Plots 1D representation or beamline elements on 2 plots: 
    ZX representation on ax_X subplot and ZY representation on ax_Y subplot
    """
    
    start_point = np.array([0,0])
    
    element_list = BL.BL_df['BL object'].values
    
    fig = plt.figure(title, figsize=(18, 12))
    
    # increase vertical space between plots
    fig.subplots_adjust(hspace=0.4)
    
    ax_X = fig.add_subplot(2,1,1)
    ax_Y = fig.add_subplot(2,1,2)
    ylim = 0.25
    
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
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='blue', alpha=0.6)
            ax_X.add_collection(p)
            
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            
            pos_x = (max(pts_x)+min(pts_x))/2
            pos_y = (max(pts_y)+min(pts_y))/2
            
            ax_X.text(pos_x, pos_y, element.name, \
                      verticalalignment='center', horizontalalignment='center')
            
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='blue', alpha=0.6)
            ax_Y.add_collection(p)
            ax_Y.text(pos_x, pos_y, element.name, \
                      verticalalignment='center', horizontalalignment='center')
        
        
        elif element.element_type == 'SM' :
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='darkblue', alpha=0.6)
            ax_X.add_collection(p) 
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            pos_x = (max(pts_x)+min(pts_x))/2
            pos_y = (max(pts_y)+min(pts_y))/2
            
            ax_X.text(pos_x, pos_y, element.name, rotation = 90,
                      verticalalignment='center', horizontalalignment='center')
            
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='darkblue', alpha=0.6)
            ax_Y.add_collection(p)
            ax_Y.text(pos_x, pos_y, element.name, rotation = 90,
                      verticalalignment='center', horizontalalignment='center')
        
        
        elif element.element_type == 'quad' :
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZX')
            p = PatchCollection(patches, color='red', alpha=0.6)
            ax_X.add_collection(p) 
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            pos_x = (max(pts_x)+min(pts_x))/2
            pos_y = (max(pts_y)+min(pts_y))/2
            
            ax_X.text(pos_x, pos_y, element.name, rotation = 90,
                      verticalalignment='center', horizontalalignment='center')
            
            [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZY')
            p = PatchCollection(patches, color='red', alpha=0.6)
            ax_Y.add_collection(p)
            ax_Y.text(pos_x, pos_y, element.name, rotation = 90,
                      verticalalignment='center', horizontalalignment='center')
            
        elif element.element_type == 'slit' :
            
            if element.orientation == 'X':
                [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZX')
                p = PatchCollection(patches, color='dimgray', alpha=0.8)
                ax_X.add_collection(p)
                # get coordinates in order to place the label
                pts_x = list(zip(*patches[0].get_xy()))[0]
                pts_y = list(zip(*patches[0].get_xy()))[1]
                pos_x = (max(pts_x)+min(pts_x))/2
                ax_X.text(pos_x, ylim*0.8, element.name, \
                          verticalalignment='center', horizontalalignment='center', \
                          rotation = 90)
                
            elif element.orientation == 'Y':
                [patches, alpha] = magnet_patches(element, start_point=start_point, orientation='ZY')
                p = PatchCollection(patches, color='dimgray', alpha=0.8)
                ax_Y.add_collection(p)
                # get coordinates in order to place the label
                pts_x = list(zip(*patches[0].get_xy()))[0]
                pts_y = list(zip(*patches[0].get_xy()))[1]
                pos_x = (max(pts_x)+min(pts_x))/2
                ax_Y.text(pos_x, ylim*0.8, element.name, \
                          verticalalignment='center', horizontalalignment='center', \
                          rotation = 90)
           
        elif element.element_type == 'BPM' :
            w = element.length
            h = 0.15
            BPM_shape = Rectangle(start_point-np.array([0,h/2]), width=w, height=h)
            p = PatchCollection([BPM_shape], color='lightblue', alpha=0.8)
            
            ax_X.add_collection(p)
            ax_X.text(start_point[0] + w/2, ylim*0.8, element.name, \
                      verticalalignment='center', horizontalalignment='center', \
                      rotation = 90)
            
            p = PatchCollection([BPM_shape], color='lightblue', alpha=0.8)
            ax_Y.add_collection(p)
            ax_Y.text(start_point[0] + w/2, ylim*0.8, element.name, \
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










def BL_geometry(BL : Beamline, refp=0):
    """ 
    Plot of the gantry structure  
    """
    
    # initial positions
    start_point = np.array([0,0])
    rot_mat = [[0,0],[0,0]]
    rot_angle = 0
    max_x = 0 # used for framing
    max_y = 0
    
    # get all elements from beamline
    element_list = BL.BL_df['BL object'].values
    
    # initiate figure
    fig = plt.figure('Beamline 2D',figsize=(18, 12))
    ax_X = fig.add_subplot(1,1,1)
    
    
    for element in element_list:
        
        if element.element_type == 'drift' and element.name!="" :
            pt = start_point + np.matmul(rot_mat, [element.length , 0])
            ax_X.plot(pt[0], pt[1] , 'o', color='black')
            ax_X.text(pt[0], pt[1]+0.1, element.name, 
                      verticalalignment='center', horizontalalignment='center')
        elif element.element_type == 'dipole' :
            # get element patches and represent it on the figure
            start_angle=rot_angle
            [patches, rot_angle] = magnet_patches(element, start_point=start_point, start_angle=start_angle, orientation='ZX', straigth_plot = False, refp=refp)
            
            p = PatchCollection(patches, color='blue', alpha=0.6)
            ax_X.add_collection(p)
            
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            # get angle at the middle of the magnet to put the label at the good position
            half_angle = (start_angle + rot_angle)/2
            text_rot_mat = [[cos(half_angle),-sin(half_angle)],[sin(half_angle),cos(half_angle)]]
            text_offset = np.matmul(text_rot_mat, [0.,0.3]) # vertical offset for text, rotated if needed
            ax_X.text(np.average(pts_x) + text_offset[0], np.average(pts_y)+text_offset[1], element.name, 
                      verticalalignment='center', horizontalalignment='center')
            
        elif element.element_type == 'SM' :
            # get element patches and represent it on the figure
            [patches, rot_angle] = magnet_patches(element, start_point=start_point, start_angle=rot_angle, orientation='ZX', straigth_plot = False)
            p = PatchCollection(patches, color='darkblue', alpha=0.6)
            ax_X.add_collection(p) 
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            text_offset = np.matmul(rot_mat, [0,0.3]) # vertical offset for text, rotated if needed
            ax_X.text(np.average(pts_x) + text_offset[0], np.average(pts_y)+text_offset[1], element.name, 
                      verticalalignment='center', horizontalalignment='center', 
                      rotation = abs(90 + rot_angle*180/3.14)%180) 
            
        
        elif element.element_type == 'quad' :
            # get element patches and represent it on the figure
            [patches, rot_angle] = magnet_patches(element, start_point=start_point, start_angle=rot_angle, orientation='ZX', straigth_plot = False)
            p = PatchCollection(patches, color='red', alpha=0.6)
            ax_X.add_collection(p) 
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            text_offset = np.matmul(rot_mat, [0,0.3]) # vertical offset for text, rotated if needed
            ax_X.text(np.average(pts_x) + text_offset[0], np.average(pts_y)+text_offset[1], element.name, 
                      verticalalignment='center', horizontalalignment='center', 
                      rotation = abs(90 + rot_angle*180/3.14)%180 )
            
            
        elif element.element_type == 'slit' :
            # get element patches and represent it on the figure
            orientation = 'Z' + element.orientation # returns ZX or ZY
            [patches, rot_angle] = magnet_patches(element, start_point=start_point, start_angle=rot_angle, orientation=orientation, straigth_plot = False)
            p = PatchCollection(patches, color='dimgray', alpha=0.8)
            ax_X.add_collection(p)   
            # get coordinates in order to place the label
            pts_x = list(zip(*patches[0].get_xy()))[0]
            pts_y = list(zip(*patches[0].get_xy()))[1]
            text_offset = np.matmul(rot_mat, [0,0.3]) # vertical offset for text, rotated if needed
            ax_X.text(np.average(pts_x) + text_offset[0], np.average(pts_y)+text_offset[1], element.name, \
                      verticalalignment='center', horizontalalignment='center', 
                      rotation = abs(90 + rot_angle*180/3.14)%180  )
            
    
        # update starting point for the next element  
        if element.curvature == 0 :
            start_point = start_point + np.matmul(rot_mat, [element.length , 0])
        else:
            rho = 1/element.curvature
            magnet_angle = element.length/rho 
            
            if element.element_rot_rad == 0 or element.element_rot_rad == np.pi or element.element_rot_rad == -np.pi :
                # find start point for next element
                # need to adjust with element_rot_rad depending if the magnet is bending upwards or downwards
                # valid only if the beamline is in a plane (i.e. element_rot_rad=0, pi or -pi)
                start_point = start_point + np.matmul(rot_mat, [rho*sin(magnet_angle) , rho*(1-cos(magnet_angle))*cos(element.element_rot_rad)])
            else:
                print("element_rot_rad = ",(np.degrees(element.element_rot_rad)))
                warnings.warn("Beamline is not in plane, no 2D representation possible")
                
            
        # update rotation matrix to represent next element in the good orientation
        rot_mat = [[cos(rot_angle),-sin(rot_angle)],[sin(rot_angle),cos(rot_angle)]]
        
        
        # rescale the plot
        if start_point[0]>max_x : max_x = start_point[0]
        if start_point[1]>max_y : max_y = start_point[1]
            
    
        
    plt.title('representation of the beamline layout')    
    ax_X.set_xlabel('length [m]')
    ax_X.set_ylabel('height [m]')
    ax_X.set_xlim([-0.5,max_x+1]) 
    ax_X.set_ylim([-0.5,max_y+1]) 
    ax_X.grid(which='major')
    
    # resize figure
    fig.set_size_inches(12, 12*(max_y+1.5)/(max_x+1.5))
    
    
    
    return [fig]



def magnet_patches(element: BL_Element, orientation = 'ZX', \
                   start_point = np.array([0,0]), start_angle=0, refp=0, \
                   straigth_plot = True):
    """ 
    Returns coordinates of element used to plot the magnet
    Possible orientations for plot: ZX or ZY
    start_point = coordinates of the starting point of the magnet (beam point)
    start_angle = initial orientation of the magnet in the beamline
    ref_p = impulsion of reference particle given the fields in the magnets (used to compute the magnets curvatures if not set previously)
    straigth_plot: straight_plot = boolean to plot a 1D or 2D layout (default=1D)
    """
    
    # initial coordinates of the magnet (beam entry point)
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
    
    
    
    if element.coil_height > 0:
        h = element.coil_height + 0.02 # add 2 cm for formers/insulation/etc
    else:
        h = 2 * a + 0.04 # default value for coil representation ; add 4 cm for formers
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
    
    
    
    if element.element_type == 'dipole' and not(straigth_plot):
        N_segments = 10  # divide into segments to represent curvature with polynomial
        if element.curvature == 0 and refp>0: 
            # modify only the curvature if not already specified
            element.set_curvature(element.Bfield, refp)
        else:
            warnings.warn("refp not set: cannot estimate bending curvature --> will lead to a 1D plot")
        bend_angle = element.get_bending_angle() # angle of dipole
    else:
        N_segments = 1
        bend_angle = 0
    
    
    # create rotation matrix if element is tilted
    alpha = start_angle
    rot_mat = [[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]]
    
    
    
    # list of points that describe the edges of the magnet shape
    # pt1 = inner line of upper magnet 
    # pt2 = outer line of upper magnet
    # pt3 = inner line of lower magnet 
    # pt4 = outer line of lower magnet 
    pt1 = [np.array([x0,y0]) + np.matmul(rot_mat, [0, a])]
    pt2 = [np.array([x0,y0]) + np.matmul(rot_mat, [-tan(CCT_angle)*h, a + h])]
    pt3 = [np.array([x0,y0]) + np.matmul(rot_mat, [0, -a])]
    pt4 = [np.array([x0,y0]) + np.matmul(rot_mat, [-tan(CCT_angle)*h, -a - h])]
    
    # segment lengths
    L_inner = L/N_segments # length of inner line - straight case
    L_outer = (L + 2 * tan(CCT_angle)*h)/N_segments # length of outer line - straight case
    if element.curvature != 0 and not(straigth_plot):
        # angle of section (curved elements)
        B_angle = bend_angle / N_segments
        
        # get curvature of the magnet edges
        rho = abs(1/element.curvature)
        
        if element.element_rot_rad == 0 or element.element_rot_rad == np.pi or element.element_rot_rad == -np.pi :
            # 2 options to code the left-right bending: 
            # 1) negative value of B_angle
            # 2) use of element_rot_rad
            # need to adjust with element_rot_rad depending if the magnet is bending upwards or downwards
            # valid only if the beamline is in a plane (i.e. element_rot_rad=0, pi or -pi)
            
            up_rot = np.sign(B_angle) * cos(element.element_rot_rad)
            B_angle = B_angle * cos(element.element_rot_rad)
            bend_angle = bend_angle * cos(element.element_rot_rad)
        else:
            print("element_rot_rad = ",(np.degrees(element.element_rot_rad)))
            warnings.warn("Beamline is not in plane, no 2D representation possible")
            
            
        rho1 = abs(rho - up_rot*a)*up_rot
        rho2 = abs(rho - up_rot*(a+h))*up_rot
        rho3 = abs(rho + up_rot*a)*up_rot
        rho4 = abs(rho + up_rot*(a+h))*up_rot
        
        
        for it in range(1, N_segments+1):
            # complete list of points
            pt1 = pt1 + [pt1[-1] + np.matmul(rot_mat, [rho1*sin(B_angle), rho1*(1-cos(B_angle))])]
            pt2 = pt2 + [pt2[-1] + np.matmul(rot_mat, [rho2*sin(B_angle), rho2*(1-cos(B_angle))])]
            pt3 = pt3 + [pt3[-1] + np.matmul(rot_mat, [rho3*sin(B_angle), rho3*(1-cos(B_angle))])]
            pt4 = pt4 + [pt4[-1] + np.matmul(rot_mat, [rho4*sin(B_angle), rho4*(1-cos(B_angle))])]
            
            # update angle and rotation matrix
            alpha = alpha + B_angle
            rot_mat = [[cos(alpha),-sin(alpha)],[sin(alpha),cos(alpha)]]
        
    else:
        L1 = L_inner
        L2 = L_outer
        L3 = L_inner
        L4 = L_outer
    
        for it in range(1, N_segments+1):
            # complete list of points
            pt1 = pt1 + [pt1[-1] + np.matmul(rot_mat, [L1, 0])]
            pt2 = pt2 + [pt2[-1] + np.matmul(rot_mat, [L2, 0])]
            pt3 = pt3 + [pt3[-1] + np.matmul(rot_mat, [L3, 0])]
            pt4 = pt4 + [pt4[-1] + np.matmul(rot_mat, [L4, 0])]
            
       
    # reverse order of points in the outer lines for a convex shape (could also be done in Polygon classs by is_convex argument)
    pt2.reverse() 
    pt4.reverse() 
    upper_pts = pt1 + pt2
    lower_pts = pt3 + pt4
    
    
    # create patches for plot
    patches = []
    
    polygon1 = Polygon(upper_pts, closed = True)
    patches.append(polygon1)
    
    polygon2 = Polygon(lower_pts, closed = True)
    patches.append(polygon2)
    
    return [patches, alpha]
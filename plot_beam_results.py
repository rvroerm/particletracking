# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 17:34:01 2019

@author: rvroerm
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fits import funcGaussian, errorFuncGaussian, funcDoubleGaussian, errorFuncDoubleGaussian
from scipy.optimize import leastsq
from transfer_functions import PtoE, Brho_scaling, split_transport_file, gaussian, transport_count_lines
from scipy.stats import norm
from BL_geometry import GTR_layout_from_transport

plt.close('all')

np.set_printoptions(precision=3)

def plot_beam(input_file,beam,it_z_ISO,it_z_GTR,ref_p):
    """
    Plot beam properties
    
    Retun spot size parameters at isocenter
    
    """
    
    
    refE = PtoE(ref_p)
    
    #######################################
    # compute GTR drawing
    nb_pts_z = transport_count_lines(input_file,1) 
    
    layout = GTR_layout_from_transport(input_file,nb_pts_z,refE)
    plt.figure('Gantry layout',figsize=(6, 4))
    plt.scatter(layout[0:nb_pts_z-1,0],layout[0:nb_pts_z-1,1])
    plt.title('Gantry layout')
    plt.xlabel('length [m]')
    plt.ylabel('heigth [m]')
    
    
    #######################################
    # Particle traces
        
    fig = plt.figure('Particle traces',figsize=(9, 6))

    fig.suptitle('Particle traces',fontweight="bold")
    
    ax0 = fig.add_subplot(211) # add subplot 1 (211 = 2 rows, 1 column, first plot)
    ax1 = fig.add_subplot(212) # add subplot 2 
    
    
    ax0.plot(beam[0:it_z_ISO,0,:],beam[0:it_z_ISO,1,:]) # plot X on left plot
    plt.ylim((-0.05,0.05))
    ax0.set_xlabel('z [m]')
    ax0.set_ylabel('x [m]')
    #plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z_ISO,6,:])))
    ax0.grid(which='major')
    
    ax1.plot(beam[0:it_z_ISO,0,:],beam[0:it_z_ISO,3,:]) # plot Y on right plot
    plt.ylim((-0.05,0.05))
    ax1.set_xlabel('z [m]')
    ax1.set_ylabel('y [m]')
    ax1.grid(which='major')
    #plt.legend(np.vectorize(PtoE)(ref_p*(1+beam[it_z_ISO,6,:])))
    
    
    
    
    ################################################
    
    plt.figure('Spot size')
    
    plt.title('Spot size at isocenter')
    
    
    
    
    # filter NaN values
    X_iso_filtered = beam[it_z_ISO,1,:][~np.isnan(beam[it_z_ISO,1,:])]
    Y_iso_filtered = beam[it_z_ISO,3,:][~np.isnan(beam[it_z_ISO,3,:])]
    
    ## add filter on central values for the fit
    #maskX = np.logical_and(X_iso_filtered < 0.01,X_iso_filtered>-0.01)
    #maskY = np.logical_and(Y_iso_filtered < 0.01,Y_iso_filtered>-0.01)
    #
    #(muX, sigmaX) = norm.fit(X_iso_filtered[maskX])
    #print('mu X = ',muX,' , sigma X = ',sigmaX)
    #(muY, sigmaY) = norm.fit(Y_iso_filtered[maskY])
    #print('mu Y = ',muY,' , sigma Y = ',sigmaY)
    
    
    nb_bins = 100
    
    n, bins, patches = plt.hist(X_iso_filtered,nb_bins,range=(-0.03,0.03),alpha=0.5,label="X")
    X_spot = [bins[:-1], n]
    
    #hist_area = np.sum(n)*(max(bins)-min(bins))/len(bins)
    ## add a 'best fit' line
    #y = norm.pdf(bins, loc=muX, scale=sigmaX)*hist_area
    #plt.plot(bins, y, 'b--', linewidth=2,label='X fit: \u03BC = {:.2f} mm, \u03C3 = {:.2f} mm'.format(muX*1000,sigmaX*1000))
    
    # make a Gaussian fit
    #tplInitial contains the "first guess" of the parameters 
    param_initial=[100,0,0.003]
    param_final,success = leastsq(errorFuncGaussian,param_initial[:],args=(X_spot[:][0],X_spot[:][1]))
    (muX,sigmaX) = param_final[1:3]
    print('mu X = %0.2e mm, sigma X = %0.2e mm'%(muX*1000,sigmaX*1000))
    
    
    plt.plot(X_spot[:][0],funcGaussian(param_final,X_spot[:][0]),'b--',label='X fit: \u03BC = {:.2f} mm, \u03C3 = {:.2f} mm'.format(muX*1000,sigmaX*1000))
    
    
    # double Gaussian fit
    param_initial=[100,0,0.003,10,0,0.01]
    param_final,success = leastsq(errorFuncDoubleGaussian,param_initial[:],args=(X_spot[:][0],X_spot[:][1]))
    print('double Gaussian X: ',param_final)
    #plt.plot(X_spot[:][0],funcDoubleGaussian(param_final,X_spot[:][0]),'b:',label='double Gaussian')
    
    
    
    n, bins, patches = plt.hist(Y_iso_filtered,nb_bins,range=(-0.03,0.03),alpha=0.5,label="Y")
    Y_spot = [bins[:-1], n]
    
    ##hist_area = np.sum(n)*(max(bins)-min(bins))/len(bins)
    ## add a 'best fit' line
    #y = norm.pdf(bins, loc=muY, scale=sigmaY)*hist_area
    #plt.plot(bins, y, 'r--', linewidth=2,label='Y fit: \u03BC = {:.2f} mm, \u03C3 = {:.2f} mm'.format(muY*1000,sigmaY*1000))
    
    param_initial=[100,0,0.003]
    param_final,success = leastsq(errorFuncGaussian,param_initial[:],args=(Y_spot[:][0],Y_spot[:][1]))
    (muY,sigmaY) = param_final[1:3]
    
    print('mu Y = %0.2e mm, sigma Y = %0.2e mm'%(muY*1000,sigmaY*1000))
    plt.plot(Y_spot[:][0],funcGaussian(param_final,Y_spot[:][0]),'r--',label='Y fit: \u03BC = {:.2f} mm, \u03C3 = {:.2f} mm'.format(muY*1000,sigmaY*1000))
    
    
    # double Gaussian fit
    param_initial=[100,0,0.003,10,0,0.01]
    param_final,success = leastsq(errorFuncDoubleGaussian,param_initial[:],args=(Y_spot[:][0],Y_spot[:][1]))
    print('double Gaussian Y: ',param_final)
    #plt.plot(Y_spot[:][0],funcDoubleGaussian(param_final,Y_spot[:][0]),'r:',label='double Gaussian')
     
    
    plt.xlabel('X/Y [m]')
    plt.ylabel('counts')
    plt.legend(loc='upper right')
    
    
    ################################################
    
    
    plt.figure('Spot size vs enrgy')
    plt.scatter(beam[it_z_ISO,1,:],beam[it_z_ISO,3,:], c = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])))
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.colorbar()
    plt.xlim((-0.05,0.05))
    plt.ylim((-0.05,0.05))
    
    
    
    
    plt.figure('Energy transmission')
    plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E source")
    plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z_GTR,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_in GTR")
    plt.hist(np.vectorize(PtoE)(ref_p*(1+beam[it_z_ISO,6,:])),100,range=(refE-10,refE+10),alpha=0.3,label="E_out GTR")
    plt.title('Energy transmission')
    plt.legend(loc='upper right')
    plt.xlabel('E [MeV]')
    plt.ylabel('nb of protons (arb. units)')
    
    
    
    ################################################
    
    fig = plt.figure('Emittance',figsize=(9, 6))
    
    
    ax0 = fig.add_subplot(121) # add subplot 1 (121 = 1 row, 2 columns, first plot)
    ax1 = fig.add_subplot(122) # add subplot 2 
    
    
    ax0.scatter(beam[it_z_ISO,1,:],beam[it_z_ISO,2,:], c = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])))
    ax0.set_xlabel('x [m]')
    ax0.set_ylabel('divx [mrad]')
    
    
    
    
    im1 = ax1.scatter(beam[it_z_ISO,3,:],beam[it_z_ISO,4,:], c = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:])))
    ax1.set_xlabel('y [m]')
    ax1.set_ylabel('divy [mrad]')
    
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation='vertical')
    
    plt.subplots_adjust(wspace=0.4) # put more space between plts
    
    return [sigmaX, sigmaY]


def get_spot_size(it_z,beam):
    """
    Get beam spot size at index it_z

    """
    
    
    X_iso_filtered = beam[it_z,1,:][~np.isnan(beam[it_z,1,:])]
    Y_iso_filtered = beam[it_z,3,:][~np.isnan(beam[it_z,3,:])]
    
    # fit
    (muX, sigmaX) = norm.fit(X_iso_filtered)
    (muY, sigmaY) = norm.fit(Y_iso_filtered)
    
    
    return [sigmaX, sigmaY]
    
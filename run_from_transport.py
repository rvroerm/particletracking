# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 11:15:31 2020

@author: rvroerm
"""


import numpy as np
from transfer_functions import transport_input, transport_count_lines, EtoP, PtoE, Brho_scaling, split_transport_file
from plot_beam_results import plot_beam, get_spot_size



def run_from_transport(input_file = "C:/TRANS/for001.dat", nb_part=1000, \
                       N_segments = 10, kill_lost_particles = True, \
                       refE = 160, old_refE = 160, DeltaE=0, E_dist='uniform',  \
                       DeltaX = 10**-5, DeltaY = 10**-5, size_dist='cst', \
                       DeltaDivX = 0.05, DeltaDivY = 0.05, div_dist='uniform', \
                       gap = 0.03, k1 = 0.5, k2 = 0, \
                       paraxial_correction = False, dpzTolerance = 10**-4, \
                       plot_results = True, output_results=True):
    
    """
    Compute trajectories of protons in beamline defined in the Transport file 'input_file'
    
    Definition of inputs:
        
    input_file = Transport input file.
    nb_part = number of particles to generate
    N_segments = number of points that represent the particle trajectory in each magnet
    kill_lost_particles = boolean to determine if paticle trajectories are stopped when they hit an element in the beamline
    refE = reference energy of the particle
    old_refE = reference energy in the Transport file (in case magnetic fields need to be scaled)
    E_dist = Statistical distribution of particle energies (uniform', 'normal' or 'cst')
    DeltaX = beam extention in X plane
    DeltaY = beam extention in Y plane
    size_dist = Statistical distribution of the beam size (uniform', 'normal' or 'cst')
    DeltaDivX = divergence in X plane
    DeltaDivY = divergence in Y plane
    div_dist = Statistical distribution of the beam divergence (uniform', 'normal' or 'cst')
    gap = vertical gap of dipoles (useful if kill_lost_particles = True)
    k1, k2 = fringe field parameters (see Transport manual)
    paraxial_correction = Boolean to correct the effective field for diverging particles (increases computation time)
    dpzTolerance = tolerance on normalized magnetic field in paraxial correction
    plot_results = boolean 
    output_results = boolean
    
    
    Definition of outputs:
        
    sigmaX, sigmaY = spot size at isocenter
    eff_ESS_dEonE_1pc, eff_GTR_dEonE_1pc = ESS and GTR effieciency for a dE/E=1% (if DeltaE high enough in input)    
    """



    split_transport_file(input_file)
    
    ref_p = EtoP(refE)
    #Brho  = 1/300*sqrt(refE**2+2*938*refE)
    
    
    Brho_factor = Brho_scaling(old_refE,refE)    
    
    
    gapX = gap # case of CCT magnets
    
    
    
    ########################################
    
    nb_pts_z = transport_count_lines(input_file,N_segments) 
    
    ###############################################################################
    # initial conditions
    
    beam = np.empty(shape=[nb_pts_z,7,nb_part]) 
    
    for i in range(0,nb_part):
        
        # initialize particle properties
        
        z=0  
        it_z = 0
        
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
        
        
        p = EtoP(E)
        dponp = (p-ref_p)/ref_p
        beam[0,:,i] = np.array([z,sizeX,divX,sizeY,divY,0,dponp])
        
    ###############################################################################    
    # extraction section
    last_z = beam[it_z,0,:]    
    last_itz = it_z 
    
    
    for i in range(0,nb_part):  
        it_z = last_itz 
        
        [beam,it_z] = transport_input('transport_file_ESS.txt',beam,refE,i,N_segments,gap,k1,k2,last_z,last_itz,Brho_factor,kill_lost_particles,gap_X=gapX,paraxial_correction=paraxial_correction,dpz_tolerance=dpzTolerance)
        
    
    
    E_list_in = np.vectorize(PtoE)(ref_p*(1+beam[0,6,:]))
    dEonE_1pc = refE/100
    nb_part_in_ESS_dEonE1pc = ((E_list_in < refE+dEonE_1pc) & (E_list_in > refE-dEonE_1pc)).sum()
    nb_part_in_ESS_dEonE2pc = ((E_list_in < refE+2*dEonE_1pc) & (E_list_in > refE-2*dEonE_1pc)).sum()
    
    E_list_out_ESS = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
    E_list_out_ESS[np.isnan(E_list_out_ESS)] = 0 
    
    
    
    
    nb_part_out_ESS_dEonE1pc = ((E_list_out_ESS < refE+dEonE_1pc) & (E_list_out_ESS > refE-dEonE_1pc)).sum()
    
    
    
    
    ###############################################################################    
    # gantry section
    last_z = beam[it_z,0,:]
    last_itz = it_z 
    it_z_GTR = it_z
    
    for i in range(0,nb_part):  
        it_z = last_itz 
        [beam,it_z] = transport_input('transport_file_GTR.txt',beam,refE,i,N_segments,gap,k1,k2,last_z,last_itz,Brho_factor,kill_lost_particles,gap_X=gapX,dpz_tolerance=dpzTolerance)
    
    
        
    # make plots
    if plot_results:
        [sigmaX, sigmaY] = plot_beam(input_file,beam,it_z,it_z_GTR,ref_p)
    else:
        [sigmaX, sigmaY] = get_spot_size(it_z,beam)
    
    
    
    # output results
    if output_results:
        
        eff_ESS_dEonE_1pc = nb_part_out_ESS_dEonE1pc/max(nb_part_in_ESS_dEonE1pc,1)*100
        print('ESS efficiency within E range [ %0.2f , %0.2f ] = %0.2f %%'%(refE-dEonE_1pc, refE+dEonE_1pc, eff_ESS_dEonE_1pc))
        
        
        E_list_out_GTR = np.vectorize(PtoE)(ref_p*(1+beam[it_z,6,:]))
        E_list_out_GTR[np.isnan(E_list_out_GTR)] = 0 
        
        
        nb_part_out_GTR_dEonE1pc = ((E_list_out_GTR < refE+dEonE_1pc) & (E_list_out_GTR > refE-dEonE_1pc)).sum()
        
        eff_GTR_dEonE_1pc = nb_part_out_GTR_dEonE1pc/max(nb_part_out_ESS_dEonE1pc,1)*100
        
        print('GTR efficiency within E range [ %0.2f , %0.2f ] = %0.2f %%'%(refE-dEonE_1pc, refE+dEonE_1pc, eff_GTR_dEonE_1pc))
        print('Total efficiency within E range [ %0.2f , %0.2f ] = %0.2f %%'%(refE-dEonE_1pc, refE+dEonE_1pc, nb_part_out_GTR_dEonE1pc/max(nb_part_in_ESS_dEonE1pc,1)*100))
        
        
        nb_part_out_GTR_dEonE2pc = ((E_list_out_GTR < refE+2*dEonE_1pc) & (E_list_out_GTR > refE-2*dEonE_1pc)).sum()
        eff_tot_dEonE_2pc = nb_part_out_GTR_dEonE2pc/max(nb_part_in_ESS_dEonE2pc,1)*100
        
        print('Total efficiency within E range [ %0.2f , %0.2f ] = %0.2f %%'%(refE-2*dEonE_1pc, refE+2*dEonE_1pc, eff_tot_dEonE_2pc))
        
    
    return [sigmaX, sigmaY, eff_ESS_dEonE_1pc, eff_GTR_dEonE_1pc]
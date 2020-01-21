# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:15:26 2019

@author: rvroerm
"""

import os
import math as math
import numpy as np
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt
from scipy.optimize import curve_fit

import warnings
import time
#import matplotlib.pyplot as plt


def transport_count_lines(transport_file,N_segments = 10,print_elements=False):
    """
    counts the number of relevant lines in file
    N_segments is the number of points displayed by element
    print_elements = boolean to print lines read in file
    """
    
    it_z = 0    
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       cnt = 1
       while line:
           
           data = line.split() 
           
           if data: # string is not empty
               if  data[0][0:2] == "4." or data[0][0:2] == "5." or data[0][0:3] == "50." or data[0][0:3] == "18.":
                   it_z = it_z + N_segments
                   if print_elements: print(data)
               elif data[0][0:2] == "3." or data[0][0:2] == "2.":
                   it_z = it_z + 1
                   if print_elements: print(data)
               elif data[0][0:5] == '(slit':
                   it_z = it_z + N_segments - 1 
                   if print_elements: print(data)
                   
           line = fp.readline()
           cnt += 1
           
    if print_elements: print('----')
    return it_z+1


def split_transport_file(transport_file):
    """
    split file between different beamline sections
    saves the output files in the current folder
    """
    
    initial_beam_found = False
    gantry_found = False
    ISO_found = False
    
    with open(transport_file) as infile, open('transport_file_ESS.txt', 'w') as outfile:
        copy = False
        for line in infile:
            if line.strip().replace(" ","") == "/*Initialbeam*/":
                initial_beam_found = True
                copy = True
            elif line.strip().replace(" ","") == "/*Gantry*/":
                gantry_found = True
                copy = False
            elif copy:
                outfile.write(line)
    
    outfile.close()
    
    with open(transport_file) as infile, open('transport_file_GTR.txt', 'w') as outfile:
        copy = False
        for line in infile:
            if line.strip().replace(" ","") == "/*Gantry*/":
                copy = True
            elif line.strip().replace(" ","") == "/*ISO*/":
                ISO_found = True
                copy = False
            elif copy:
                outfile.write(line)
     
    
    if not initial_beam_found:
        warnings.warn(" /*Initialbeam*/ label not defined in Transport file, plots might not correspond to the desired location")            
    if not gantry_found:
        warnings.warn(" /*Gantry*/ label not defined in Transport file, plots might not correspond to the desired location")
    if not ISO_found:
        warnings.warn(" /*ISO*/ label not defined in Transport file, plots might not correspond to the desired location")
        
    outfile.close()



def Brho_scale_transport(transport_file,newE,new_file_path="D:/temp/",rest_mass=938):
    """
    B rho scaling from a transport file
    code will fail if no line "1.xxxx" is found before other elements (which should never happen)
    """
    
    filename = os.path.basename(transport_file)
    filename_noext = os.path.splitext(filename)[0] + "_" + str(newE) + "MeV" + ".txt"
    
    
    with open(transport_file) as infile, open(new_file_path + filename_noext, 'w') as outfile: 
       line = infile.readline()
       cnt = 1
       
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           
           if data: # string is not empty
               
               if data[0][0:2] == "1.":
                   oldE = PtoE(float(data[7].replace(";","") ))
                   scaling_ratio = Brho_scaling(oldE,newE,rest_mass)
                   data[7] = EtoP(newE)
                   outfile.write("".join(str(data)).replace(","," ").replace('[',"").replace(']',"").replace("'","") + '\n')
                   
               elif data[0][0:2] == "4.":
                   data[2] = float(data[2].replace(";","") )*scaling_ratio
                   outfile.write("".join(str(data)).replace(","," ").replace('[',"").replace(']',"").replace("'","") + '\n')
                   
               elif data[0][0:2] == "5.":
                   data[2] = float(data[2].replace(";","") )*scaling_ratio   
                   outfile.write("".join(str(data)).replace(","," ").replace('[',"").replace(']',"").replace("'","") + '\n')
                   
               elif data[0][0:2] == "18.":
                   data[2] = float(data[2].replace(";","") )*scaling_ratio
                   outfile.write("".join(str(data)).replace(","," ").replace('[',"").replace(']',"").replace("'","") + '\n')
                   
               else:
                   outfile.write(line)
           else:
               outfile.write(line)
             
           line = infile.readline()
           cnt += 1
           
           
    

def transport_input(transport_file,beam,refE,it_p,N_segments,gap,k1,k2,z,it_z,scaling_ratio = 1,kill_lost_particles=True,gap_X=1,paraxial_correction = False,dpz_tolerance=10**-6):
    """
    opens a transport file and computes beam through magnets for one particle
    scaling_ratio can be used to make a Brho scaling of the transport file
    """
    
    # default settings
    B = 0
    
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       prev_line = fp.readline()
       cnt = 1
       
       
       
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           data_prev = prev_line.split()
           
           
           
           if data: # string is not empty
               
               if data[0][0:2] == "2.":
                   if data_prev :
                       if data_prev[0] == "4.":
                           # exit pole face (entrance managed within 4.0)
                           angle = float(data[1].replace(";","") )
                           [beam[:,:,it_p],it_z] = pole_face(angle,B,gap,k1,k2,beam[:,:,it_p],refE,it_z, paraxial_correction = paraxial_correction)
                           
                           
               if data[0][0:2] == "3.":
                   # drift
                   L = float(data[1].replace(";","") )
                   z = z + L
                   
                   [beam[:,:,it_p],it_z] = drift(L,beam[:,:,it_p],refE,it_z,1)
                   
               if data[0][0:2] == "4.":
                   # sector bending
                   
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )*scaling_ratio
                   n = float(data[3].replace(";","") )
                   
                   
                   if data_prev :
                       if data_prev[0] == "2.": # entrance pole face
                           angle = float(data_prev[1].replace(";","") )
                           [beam[:,:,it_p],it_z] = pole_face(angle,B,gap,k1,k2,beam[:,:,it_p],refE,it_z, paraxial_correction = paraxial_correction)
                           
                   
                   
                   [beam[:,:,it_p],it_z] = sec_bend(L,B,n,gap,beam[:,:,it_p],refE,it_z,N_segments=N_segments, \
                   kill_lost_particles=kill_lost_particles,gap_X=gap_X,paraxial_correction = paraxial_correction, dpz_tolerance=dpz_tolerance)
                   
                   
               if data[0][0:2] == "5.": 
                   # quad
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )*scaling_ratio
                   a = float(data[3].replace(";","") )/1000 # aperture converted to meters
                   
                   #[beam[:,:,it_p],it_z] = quad(L,B,a,beam[:,:,it_p],it_z,refE,N_segments,kill_lost_particles)
                   
                   [beam[:,:,it_p],it_z] = quad(L,B,a,beam[:,:,it_p],it_z,refE,N_segments=N_segments, \
                   kill_lost_particles=kill_lost_particles,paraxial_correction=paraxial_correction,dpz_tolerance=dpz_tolerance)
                   
               
               if data[0][0:3] == "18.": 
                   # sextupole
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )*scaling_ratio
                   a = float(data[3].replace(";","") )/1000 # aperture converted to meters
                   angle = float(data[4].replace(";","") ) # not used in original transport
                   
                   [beam[:,:,it_p],it_z] = sextupole(L,B,a,angle,beam[:,:,it_p],it_z,refE,N_segments,kill_lost_particles) 
                
                
               if data[0][0:3] == '50.':
                   # linear wedge
                   L = float(data[1].replace(";","") )
                   width = float(data[2].replace(";","") )
                   mat = data[3].replace(";","") 
                   
                   [beam[:,:,it_p],it_z] = wedge_X(L,width,mat,beam[:,:,it_p],it_z,refE,N_segments)
                   
               if data[0][0:5] == '(slit':
                   # slit
                   
                   L = float(data[1].replace(";","") )
                   leftX = float(data[2].replace(";","") )
                   rightX = float(data[3].replace(";","") )
                   mat = data[4].replace(";","") 
                   
                   [beam[:,:,it_p],it_z] = slitX(L,leftX,rightX,mat,beam[:,:,it_p],it_z,refE,N_segments)
                   
                   line = fp.readline() # skip next drift
                   cnt +=1 
                   
                   
           prev_line = line     
           line = fp.readline()
           cnt += 1
    
    
    return [beam, it_z]




def drift(L,beam,refE,it_z,N_segments=1):
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    gamma = PtoGamma(p,938)
    
    
    
    drift_mat = np.array([[1,L_segment,0,0,0,0],[0,1,0,0,0,0],[0,0,1,L_segment,0,0],[0,0,0,1,0,0],[0,0,0,0,1,L_segment/gamma**2],[0,0,0,0,0,1]])
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + L_segment
        beam[it_z,1:7] = np.matmul(drift_mat,np.transpose(beam[it_z-1,1:7]))
        
        L_eff = abs(L_segment/cos(sqrt(beam[it_z,2]**2+beam[it_z,4]**2))) # effective length travelled by particle (negliect change in direction)
        beam[it_z,5] = (L_segment - L_eff) + beam[it_z,5] #higher order correction related to beam angle
        
    return [beam, it_z]




def drift_core(L,beam,refE):
    
    ref_p = EtoP(refE)
    p = ref_p + beam[5]*ref_p
    gamma = PtoGamma(p,938)
    
    drift_mat = np.array([[1,L,0,0,0,0],[0,1,0,0,0,0],[0,0,1,L,0,0],[0,0,0,1,0,0],[0,0,0,0,1,L/gamma**2],[0,0,0,0,0,1]])
    beam = np.matmul(drift_mat,np.transpose(beam))
    
    return beam

def drift_matrix(L,p,rest_mass=938):
    gamma = PtoGamma(p,938)
    drift_mat = np.array([[1,L,0,0,0,0],[0,1,0,0,0,0],[0,0,1,L,0,0],[0,0,0,1,0,0],[0,0,0,0,1,L/gamma**2],[0,0,0,0,0,1]])
    return drift_mat    

def RF_cavity(delta_E,cav_length,freq,phi,beam,refE,it_z,N_segments = 10):
    # simplified version, neglects defocusing effects
    
    
    omega = 2*math.pi*freq
    ref_p = EtoP(refE)
    
    for i in range(0,N_segments):
        
        p = ref_p + beam[it_z,6]*ref_p
        speed_ref = PtoV(ref_p,938)
        
        E_i = PtoE(p)
        
        z_eff = beam[it_z,0] + beam[it_z,5]
        
        t = z_eff/speed_ref # divide by speed_ref because speed difference is included in z_eff
        
        E_f = E_i + delta_E/N_segments*cos(omega*t+phi) 
        
        p_f = EtoP(E_f)
        
        [beam,it_z] = drift(cav_length/N_segments,beam,refE,it_z,1)
        
        beam[it_z,6] = (p_f - ref_p)/ref_p
        
    
    return [beam, it_z]

def RF_cavity_old(delta_E,freq,phi,beam,refE,it_z):
    # simplified version, neglects defocusing effects
    
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    speed_ref = PtoV(ref_p,938)
    
    
    E_i = PtoE(p)
    
    z_eff = beam[it_z,0] + beam[it_z,5]
    
    t = z_eff/speed_ref
    
    omega = 2*math.pi*freq
    
    E_f = E_i + delta_E*cos(omega*t+phi) 
    
    p_f = EtoP(E_f)
    
    #print('E_i = ',E_i,', E_f = ',E_f,', t = ',t)
    #print(beam[it_z,:])
    
    it_z += 1
    beam[it_z,0:5] = beam[it_z-1,0:5]
    beam[it_z,6] = (p_f - ref_p)/ref_p
    
    
    return [beam, it_z]
    






def quad(L,B,a,beam,it_z,refE,N_segments = 10,kill_lost_particles=True,rest_mass=938, \
          paraxial_correction = False, dpz_tolerance=10**-6):
    """
    L =length, B=field, a=aperture
    kill_lost_particles=True will kill particles outside the good field region
    paraxial_correction = True will scale the fields according to take into account the fact that pz != p
    the relative tolerace is dpz_tolerance
    """
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    
    
    for i in range(0,N_segments):
        
        if np.isnan(beam[it_z,1]) or (kill_lost_particles and (sqrt(beam[it_z,1]**2 + beam[it_z,3]**2) >= a)):
            # lost particle
            it_z = it_z + 1
            beam[it_z,0] = beam[it_z-1,0] + L_segment
            beam[it_z,1:7] = np.nan            
        else:
            if paraxial_correction :
                corr_fact = sqrt(beam[it_z,2]**2+beam[it_z,4]**2+1)
                B_corrected = B / corr_fact
            else:
                corr_fact = 1
                B_corrected = B
            
            quad_mat = quad_matrix(L_segment,B_corrected,a,p,rest_mass)
            
            new_beam = np.matmul(quad_mat,np.transpose(beam[it_z,1:7]))
            
            
            
            if abs(sqrt(new_beam[1]**2+new_beam[3]**2+1) - corr_fact) > dpz_tolerance and paraxial_correction :
                # need more intermediate points
                
                intermediate_segments = int(abs(sqrt(new_beam[1]**2+new_beam[3]**2+1) - corr_fact) / dpz_tolerance ) + 1
                
                
                new_beam = beam[it_z,1:7]
                for j in range(0,intermediate_segments):
                    corr_fact = sqrt(new_beam[1]**2+new_beam[3]**2+1)
                    B_corrected = B / corr_fact
                    quad_mat = quad_matrix(L_segment/intermediate_segments,B_corrected,a,p,rest_mass)
                    
                    new_beam = np.matmul(quad_mat,np.transpose(new_beam))
                
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                beam[it_z,1:7] = new_beam
            else:
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                beam[it_z,1:7] = np.matmul(quad_mat,np.transpose(beam[it_z-1,1:7]))
                
                L_eff = abs(L_segment/cos(sqrt(beam[it_z,2]**2+beam[it_z,4]**2))) # effective length travelled by particle (neglect change in direction)
                beam[it_z,5] = (L_segment - L_eff) + beam[it_z,5] #higher order correction related to beam angle

            
    
    p = ref_p + beam[it_z,6]*ref_p
    
    return [beam, it_z]
    
def quad_matrix(L,B,a,p,rest_mass=938):
    """
    returns the transfer matrix for a particle of momentum p
    L =length, B=field, a=aperture
    """
    
    gamma = PtoGamma(p,rest_mass)
    
    Brho  = PtoBrho(p)
    
    #Brho  = 1/300*sqrt(refE**2+2*938*refE)
    k_q  = sqrt(abs(B)/a*(1/Brho))
    
    
    if B>0: 
        # X focusing
        quad_mat = np.array([[cos(k_q*L),1/k_q*sin(k_q*L),0,0,0,0],
                              [-k_q*sin(k_q*L),cos(k_q*L),0,0,0,0],
                            [0,0,cosh(k_q*L),1/k_q*sinh(k_q*L),0,0],
                            [0,0,k_q*sinh(k_q*L),cosh(k_q*L),0,0],
                            [0,0,0,0,1,L/gamma**2],
                            [0,0,0,0,0,1]])
    else:
        # Y focusing
        quad_mat = np.array([[cosh(k_q*L),1/k_q*sinh(k_q*L),0,0,0,0],
                              [k_q*sinh(k_q*L),cosh(k_q*L),0,0,0,0],
                            [0,0,cos(k_q*L),1/k_q*sin(k_q*L),0,0],
                            [0,0,-k_q*sin(k_q*L),cos(k_q*L),0,0],
                            [0,0,0,0,1,L/gamma**2],
                            [0,0,0,0,0,1]])
    
    return quad_mat

def sextupole(L,B,a,alpha,beam,it_z,refE,N_segments = 10):
    """
    sextupole: L =length, B=field, a=aperture, alpha = angle
    """
    
    alpha = math.radians(alpha)
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    E = PtoE(p)
    
    Brho  = 1/300*sqrt(E**2+2*938*E)
    #Brho  = 1/300*sqrt(refE**2+2*938*refE)
    k_s  = sqrt(abs(B)/a**2*(1/Brho))
    
    R11 = 1
    R12 = L_segment
    T111 = -1/2 * k_s**2 * L_segment**2
    T112 = -1/3 * k_s**2 * L_segment**3
    T122 = -1/12 * k_s**2 * L_segment**4
    T133 = - T111
    T134 = - T112
    T144 = - T122
    
    R22 = 1
    T211 = - k_s**2 * L_segment
    T212 = - k_s**2 * L_segment**2
    T222 = -1/3 * k_s**2 * L_segment**3
    T233 = - T211 
    T234 = - T212 
    T244 = - T222 
    
    R33 = 1
    R34 = L_segment
    T313 = k_s**2 * L_segment**2
    T314 = 1/3 * k_s**2 * L_segment**3
    T323 = 1/3 * k_s**2 * L_segment**3
    T324 = 1/6 * k_s**2 * L_segment**4
    
    R44 = 1
    T413 = 2 * k_s**2 * L_segment
    T414 = k_s**2 * L_segment**2
    T423 = T414
    T424 = 2/3 * k_s**2 * L_segment**3
    
    beam_rot_out = np.empty((6,1)) 
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + L_segment
        
        beam_rot_in = rotation(beam[it_z-1,1:7],alpha)
        
        beam_rot_out[0] = R11*beam_rot_in[0] + R12*beam_rot_in[1] + T111*beam_rot_in[0]**2 + T112*beam_rot_in[0]*beam_rot_in[1] + T122*beam_rot_in[1]**2 + T133*beam_rot_in[2]**2 + T134*beam_rot_in[2]*beam_rot_in[3] + T144*beam_rot_in[3]**2
        beam_rot_out[1] = R22*beam_rot_in[1] + T211*beam_rot_in[0]**2 + T212*beam_rot_in[0]*beam_rot_in[1] + T222*beam_rot_in[1]**2 + T233*beam_rot_in[2]**2 + T234*beam_rot_in[2]*beam_rot_in[3] + T244*beam_rot_in[3]**2
        beam_rot_out[2] = R33*beam_rot_in[2] + R34*beam_rot_in[3] + T313*beam_rot_in[0]*beam_rot_in[2] + T314*beam_rot_in[0]*beam_rot_in[3] + T323*beam_rot_in[1]*beam_rot_in[2] + T324*beam_rot_in[1]*beam_rot_in[3]
        beam_rot_out[3] = R44*beam_rot_in[3] + T413*beam_rot_in[0]*beam_rot_in[2] + T414*beam_rot_in[0]*beam_rot_in[3] + T423*beam_rot_in[1]*beam_rot_in[2] + T424*beam_rot_in[1]*beam_rot_in[3]
        beam_rot_out[4] = beam_rot_in[4]
        beam_rot_out[5] = beam_rot_in[5]
        
        
        beam_rot_out = rotation(beam_rot_out,-alpha)
        beam[it_z,1:7] = np.transpose(beam_rot_out)
    
    return [beam, it_z]


def sec_bend_matrix(L,B,n,p,rest_mass=938):
    
    Brho  = PtoBrho(p)
    h = B/Brho
    
    
    if n>0 and n<1:
        kx  = sqrt((1-n)*h**2)
        ky  = sqrt(n*h**2) 
        bend_mat = np.array([[cos(kx*L),1/kx*sin(kx*L),0,0,0,h/kx**2*(1-cos(kx*L))],
                              [-kx*sin(kx*L),cos(kx*L),0,0,0,h/kx*sin(kx*L)],
                              [0,0,cos(ky*L),1/ky*sin(ky*L),0,0],
                              [0,0,-ky*sin(ky*L),cos(ky*L),0,0],
                              [-h/kx*sin(kx*L),-h/kx**2*(1-cos(kx*L)),0,0,1,-h**2/kx**3*(kx*L-sin(kx*L))],
                              [0,0,0,0,0,1]])
    elif n>1:
        kx  = sqrt(-(1-n)*h**2) #complex value --> change sign and forumulas
        ky  = sqrt(n*h**2) 
        bend_mat = np.array([[cosh(kx*L),1/kx*sinh(kx*L),0,0,0,-h/kx**2*(1-cosh(kx*L))],
                              [kx*sinh(kx*L),cosh(kx*L),0,0,0,h/kx*sinh(kx*L)],
                              [0,0,cos(ky*L),1/ky*sin(ky*L),0,0],
                              [0,0,-ky*sin(ky*L),cos(ky*L),0,0],
                              [-h/kx*sinh(kx*L),h/kx**2*(1-cosh(kx*L)),0,0,1,h**2/kx**3*(kx*L-sinh(kx*L))],
                              [0,0,0,0,0,1]])
    elif n<0:
        kx  = sqrt((1-n)*h**2)
        ky  = sqrt(-n*h**2) #complex value --> change sign and forumulas
        bend_mat = np.array([[cos(kx*L),1/kx*sin(kx*L),0,0,0,h/kx**2*(1-cos(kx*L))],
                              [-kx*sin(kx*L),cos(kx*L),0,0,0,h/kx*sin(kx*L)],
                              [0,0,cosh(ky*L),1/ky*sinh(ky*L),0,0],
                              [0,0,ky*sinh(ky*L),cosh(ky*L),0,0],
                              [-h/kx*sin(kx*L),-h/kx**2*(1-cos(kx*L)),0,0,1,-h**2/kx**3*(kx*L-sin(kx*L))],
                              [0,0,0,0,0,1]])
    elif n==1: #kx=0
        ky  = sqrt(n*h**2) 
        bend_mat = np.array([[1,L,0,0,0,0],
                              [0,1,0,0,0,0],
                              [0,0,cos(ky*L),1/ky*sin(ky*L),0,0],
                              [0,0,-ky*sin(ky*L),cos(ky*L),0,0],
                              [-h*L,-h*L**2/2,0,0,1,-h**2*L**3/6],
                              [0,0,0,0,0,1]])
    elif n==0:
        kx  = sqrt((1-n)*h**2) #ky=0
        bend_mat = np.array([[cos(kx*L),1/kx*sin(kx*L),0,0,0,h/kx**2*(1-cos(kx*L))],
                              [-kx*sin(kx*L),cos(kx*L),0,0,0,h/kx*sin(kx*L)],
                              [0,0,1,L,0,0],
                              [0,0,0,1,0,0],
                              [-h/kx*sin(kx*L),-h/kx**2*(1-cos(kx*L)),0,0,1,-h**2/kx**3*(kx*L-sin(kx*L))],
                              [0,0,0,0,0,1]])
    
    
    
    
    return bend_mat





def sec_bend(L,B,n,gap,beam,refE,it_z,N_segments = 10,kill_lost_particles=True,rest_mass=938,gap_X=1, \
          paraxial_correction = False, dpz_tolerance=10**-6):
    """
    L =length, B=field, a=aperture
    kill_lost_particles=True will kill particles outside the good field region
    paraxial_correction = True will scale the fields according to take into account the fact that pz != p
    the relative tolerace is dpz_tolerance
    """
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    
    
    for i in range(0,N_segments):
        
        if np.isnan(beam[it_z,1]) or (kill_lost_particles and (abs(beam[it_z,1]) >= gap_X or abs(beam[it_z,3]) >= gap)):
            # lost particle
            it_z = it_z + 1
            beam[it_z,0] = beam[it_z-1,0] + L_segment
            beam[it_z,1:7] = np.nan            
        else:
            if paraxial_correction :
                corr_fact = sqrt(beam[it_z,2]**2+beam[it_z,4]**2+1)
                B_corrected = B / corr_fact
            else:
                corr_fact = 1
                B_corrected = B
            
            bend_mat = sec_bend_matrix(L_segment,B_corrected,n,p,rest_mass)
            
            new_beam = np.matmul(bend_mat,np.transpose(beam[it_z,1:7]))
            
            
            
            if abs(sqrt(new_beam[1]**2+new_beam[3]**2+1) - corr_fact) > dpz_tolerance and paraxial_correction :
                # need more intermediate points
                
                intermediate_segments = int(abs(sqrt(new_beam[1]**2+new_beam[3]**2+1) - corr_fact) / dpz_tolerance ) + 1
                
                new_beam = beam[it_z,1:7]
                for j in range(0,intermediate_segments):
                    corr_fact = sqrt(new_beam[1]**2+new_beam[3]**2+1)
                    B_corrected = B / corr_fact
                    
                    bend_mat = sec_bend_matrix(L_segment/intermediate_segments,B_corrected,n,p,rest_mass)
                    
                    new_beam = np.matmul(bend_mat,np.transpose(new_beam))
                
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                beam[it_z,1:7] = new_beam
            else:
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                beam[it_z,1:7] = np.matmul(bend_mat,np.transpose(beam[it_z-1,1:7]))
                
            
    
    p = ref_p + beam[it_z,6]*ref_p
    
    return [beam, it_z]


def pole_face(angle,B,gap,k1,k2,beam,refE,it_z,paraxial_correction = False):
    """
    k1 and k2 related to frindge field properties; default k1=0.5, k2=0
    see SLAC transport manual
    """
    
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    Brho  = PtoBrho(p)
    
    if paraxial_correction :
        corr_fact = sqrt(beam[it_z,2]**2+beam[it_z,4]**2+1)
        B_corrected = B / corr_fact
    else:
        B_corrected = B
    
    
    radius = Brho/B_corrected
    
    beta = math.radians(angle)
    
    
    psi = k1*(gap/radius)*((1+sin(beta)**2)/cos(beta))*(1-k1*k2*gap/radius*tan(beta))    
    
    
    pole_mat = np.array([[1,0,0,0,0,0],
                    [tan(beta)/radius,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,-tan(beta-psi)/radius,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    
    it_z = it_z + 1
    beam[it_z,0] = beam[it_z-1,0]
    beam[it_z,1:7] = np.matmul(pole_mat,np.transpose(beam[it_z-1,1:7]))
    
    
    
    
    return [beam, it_z]


def pole_face_mat(angle,B,gap,p,k1=0.5,k2=0,paraxial_correction = False):
    """
    k1 and k2 related to frindge field properties; default k1=0.5, k2=0
    see SLAC transport manual
    """
    
    Brho  = PtoBrho(p)
    
    radius = Brho/B
    
    beta = math.radians(angle)
    
    
    psi = k1*(gap/radius)*((1+sin(beta)**2)/cos(beta))*(1-k1*k2*gap/radius*tan(beta))    
    
    
    pole_mat = np.array([[1,0,0,0,0,0],
                    [tan(beta)/radius,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,-tan(beta-psi)/radius,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    return pole_mat

def bending(L,B,n,pole_in,pole_out,gap,k1,k2,beam,refE,it_z,N_segments = 10,kill_lost_particles=True,gap_X=1,paraxial_correction = False, dpz_tolerance=10**-6):
    """
    kill_lost_particles will kill particles that are beyond the gap opening in Y
    horizontal gap is used to kill particles that are too far; default = 1 m
    paraxial_correction = True will scale the fields according to take into account the fact that pz != p
    the relative tolerace is dpz_tolerance
    """
    
    [beam, it_z] = pole_face(pole_in,B,gap,k1,k2,beam,refE,it_z,paraxial_correction = paraxial_correction)
    
    #[beam, it_z] = sec_bend(L,B,n,gap,beam,refE,it_z,N_segments,kill_lost_particles,gap_X)
    
    [beam, it_z] = sec_bend(L,B,n,gap,beam,refE,it_z,N_segments = N_segments,kill_lost_particles=kill_lost_particles,gap_X=gap_X, \
          paraxial_correction = paraxial_correction, dpz_tolerance=dpz_tolerance)
    
    [beam, it_z]= pole_face(pole_out,B,gap,k1,k2,beam,refE,it_z,paraxial_correction = paraxial_correction)

    return [beam, it_z]





def solenoid(L,B,a,beam,it_z,refE,N_segments=20):
    
    
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    gamma = PtoGamma(p,938)
    #speed = PtoV(p,938)
    
    
    
    
    E = PtoE(p)
    Brho  = 1/300*sqrt(E**2+2*938*E)
    
    K = 1/2 * B/Brho
    L_segment = L/N_segments
    
    
#    if Brho/B < L*10 :
#        # need more segments
#        N_segments = int(L/(Brho/B))*10
#        L_segment = L/N_segments
#    else:
#        N_segments = 1
#        L_segment = L
    

    
    entrance_mat = np.array([[1,0,0,0,0,0],
                    [0,1,K,0,0,0],
                    [0,0,1,0,0,0],
                    [-K,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    
    exit_mat = np.array([[1,0,0,0,0,0],
                    [0,1,-K,0,0,0],
                    [0,0,1,0,0,0],
                    [K,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    C = cos(K*L_segment)
    S = sin(K*L_segment)
    
    segment_mat = np.array([[1,S*C/K,0,S*S/K,0,0],
                    [0,2*C*C-1,0,2*S*C,0,0],
                    [0,-S*S/K,1,S*C/K,0,0],
                    [0,-2*S*C,0,2*C*C-1,0,0],
                    [0,0,0,0,1,L_segment/gamma**2],
                    [0,0,0,0,0,1]])
    
    it_z = it_z + 1
    beam[it_z,0] = beam[it_z-1,0]
    beam[it_z,1:7] = np.matmul(entrance_mat,np.transpose(beam[it_z-1,1:7]))
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + L_segment
        beam[it_z,1:7] = np.matmul(segment_mat,np.transpose(beam[it_z-1,1:7]))
        L_eff = abs(L_segment/cos(sqrt(beam[it_z,2]**2+beam[it_z,4]**2))) # effective length travelled by particle (negliect change in direction)
        beam[it_z,5] = (L_segment - L_eff) + beam[it_z,5] #higher order correction related to beam angle
    
    it_z = it_z + 1
    beam[it_z,0] = beam[it_z-1,0]
    beam[it_z,1:7] = np.matmul(exit_mat,np.transpose(beam[it_z-1,1:7]))
        
    
    
    return [beam, it_z]


    


def wedge_X(thickness_max,width,mat,beam,it_z,refE,N_segments=10):
    
    x = beam[it_z,1]
    
    
    a = thickness_max/width
    b = thickness_max/2
    
    thickness = max(a*x+b,0)
    
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + thickness/N_segments
        [beam[it_z,1:7],sigma] = Highland(thickness/N_segments,mat,beam[it_z-1,1:7],refE)
    
    return [beam, it_z]





def wedge_circular(L_min,a_min,L_max,a_max,mat,N_steps,beam,it_z,refE):
    # wedge made of stack of N_steps collimators, with a min thickness L_min for an aperture a_min and a max thickness L_max at an aperture a_max
    
    #entrance edge
    if abs(L_max-L_min) > 0 :
        for i in range(0,N_steps):
            a = (a_max*(N_steps-i) + a_min*i)/N_steps
            L = (L_max - L_min)/(N_steps)/2 
            #print("i = ",i,"a = ",a,", L = ",L)
            [beam, it_z] = collimator(L,a,mat,beam,it_z,1,refE) 
    
    #central part
    if L_min > 0 :
        [beam, it_z] = collimator(L_min,a,mat,beam,it_z,1,refE) 
    
    #exit edge
    if abs(L_max-L_min) > 0 :
        for i in range(N_steps-1,-1,-1):
            a = (a_max*(N_steps-i) + a_min*i)/N_steps
            L = (L_max - L_min)/(N_steps)/2 
            #print("i = ",i,"a = ",a,", L = ",L)
            [beam, it_z] = collimator(L,a,mat,beam,it_z,1,refE) 
        
    
    return beam  



def slit(orientation,L_m,left,right,mat,beam,it_z,refE,N_segments=10):
    """
    slit with 2 blades (left and right)
    """
    
    L_segment = L_m/N_segments
    
    if orientation == 'X':
        for i in range(0,N_segments):
            if (beam[it_z,1] <= left or beam[it_z,1] >= right):
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                [beam[it_z,1:7],sigma] = Highland(L_segment,mat,beam[it_z-1,1:7],refE)
            else:
                [beam,it_z] = drift(L_segment,beam,refE,it_z,1)
    elif orientation == 'Y':
        for i in range(0,N_segments):
            if (beam[it_z,3] <= left or beam[it_z,3] >= right):
                it_z = it_z + 1
                beam[it_z,0] = beam[it_z-1,0] + L_segment
                [beam[it_z,1:7],sigma] = Highland(L_segment,mat,beam[it_z-1,1:7],refE)
            else:
                [beam,it_z] = drift(L_segment,beam,refE,it_z,1)
                
    
    return [beam, it_z]


def slitX(L_m,leftX,rightX,mat,beam,it_z,refE,N_segments=10):
    """
    slit along X axis, with 2 blades (left and right)
    """
    
    L_segment = L_m/N_segments
    
    for i in range(0,N_segments):
        if (beam[it_z,1] <= leftX or beam[it_z,1] >= rightX):
            it_z = it_z + 1
            beam[it_z,0] = beam[it_z-1,0] + L_segment
            [beam[it_z,1:7],sigma] = Highland(L_segment,mat,beam[it_z-1,1:7],refE)
            
        else:
            [beam,it_z] = drift(L_segment,beam,refE,it_z,1)
    
    return [beam, it_z]

def collimator(L_m,a,mat,beam,it_z,refE,N_segments):
    """
    L_m=collimator thickness, a=aperture radius, mat=material
    """
    
    L_segment = L_m/N_segments
    
    for i in range(0,N_segments):
        if sqrt(beam[it_z,1]**2 + beam[it_z,3]**2) >= a:
            it_z = it_z + 1
            beam[it_z,0] = beam[it_z-1,0] + L_segment
            [beam[it_z,1:7],sigma] = Highland(L_segment,mat,beam[it_z-1,1:7],refE)
        else:
            [beam,it_z] = drift(L_segment,beam,refE,it_z,1) # drift only
    
    
#    if sqrt(beam[it_z,1]**2 + beam[it_z,3]**2) >= a:
#        L_segment = L_m/N_segments
#        for i in range(0,N_segments):
#            it_z = it_z + 1
#            beam[it_z,0] = beam[it_z-1,0] + L_segment
#            [beam[it_z,1:7],sigma] = Highland(L_segment,mat,beam[it_z-1,1:7])
#        
#    else:
#        [beam,it_z] = drift2(L_m,beam,it_z,N_segments) # drift only
    
    return [beam, it_z]

def rotation(beam,angle):
    
    rot_mat = np.array([[cos(angle),0,sin(angle),0,0,0],
                    [0,cos(angle),0,sin(angle),0,0],
                    [-sin(angle),0,cos(angle),0,0,0],
                    [0,-sin(angle),0,cos(angle),0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    beam = np.matmul(rot_mat,beam)
    
    return beam
    
def Highland(L_m,mat,beam_in,refE):
    """computes scattering angle"""

    ref_p = EtoP(refE)
    p = ref_p + beam_in[5]*ref_p
    E = PtoE(p)
    
    if math.isnan(E):
        beam_out = beam_in
        sigma = float('nan') 
    elif L_m <= 0 :
        beam_out = beam_in
        sigma = 0
    else:
        p = sqrt(2*938*E+E**2)
        gamma = (938+E)/938
        beta = sqrt(1-1/gamma**2)
        
        if mat == 'water':
            Lrad = 36.08
            density = 1
            cWET = 1
        elif mat == 'tantalum':
            Lrad = 6.82
            density = 16.69
            cWET = 0.5057
        elif mat == 'lexan':
            Lrad = 41.72
            density = 1.2
            cWET = 0.95237
        elif mat == 'carbon':
            Lrad = 42.7
            density = 2
            cWET = 0.894
        elif mat == 'iron':
            Lrad = 13.84
            density = 7.874
            cWET = 0.4349
        else:
            raise Exception('material not found')
        
        
        
        # modify energy
        p_range = EtoR(E)
        WET = L_m*100*cWET*density/cos(sqrt(beam_in[1]**2+beam_in[3]**2))
        
        new_p_range = p_range - WET
        if new_p_range>0:
            # assume drift is independent of the scattering (i.e. impact on angle only)
            beam_out = drift_core(L_m,beam_in,refE) 
            
            
            L_gcm2 = abs(L_m*100*density/cos(sqrt(beam_in[1]**2+beam_in[3]**2)))
            
            sigma = 14.1/p/beta*sqrt(L_gcm2/Lrad)*(1+1/9*log10(L_gcm2/Lrad))
            angleX = np.random.normal(0,sigma) 
            angleY = np.random.normal(0,sigma)
            beam_out[1] = beam_in[1] + angleX
            beam_out[3] = beam_in[3] + angleY
            
            
            
            new_E = RtoE(new_p_range)
            new_p = EtoP(new_E)
            
            # make the same calculations for the initial energy to avoid rounding errors
            old_E = RtoE(p_range)
            old_p = EtoP(old_E)
            
            beam_out[5] = beam_out[5] + (new_p-old_p)/ref_p
            
        else:
            beam_out = np.empty(6)
            beam_out[:] = np.nan
            sigma = 0
             
    

    return [beam_out, sigma]


def PtoGamma(p,rest_mass=938):
    """impulsion (in Mev/c) to gamma factor"""
    
    beta = PtoBeta(p,rest_mass)
    gamma = BetatoGamma(beta)
    return gamma


def PtoV(p,rest_mass=938):
    """impulsion (in Mev/c) to speed (in m/s)"""
    
    light_speed = 2.998*10**8
    speed = PtoBeta(p,rest_mass)*light_speed
    
    return speed

def PtoBeta(p,rest_mass=938):
    """ impulsion (in Mev/c) to speed (in 1/c)"""
    beta = p/sqrt(rest_mass**2 + p**2)
    return beta


def VtoGamma(v):
    return 1/sqrt(1-v**2/(2.998*10**8)**2)   

def GammatoV(gamma):
    return 2.998*10**8 * sqrt(1 - 1/gamma**2)

def BetatoGamma(beta):
    return 1/sqrt(1-beta**2)   

def GammatoBeta(gamma):
    return sqrt(1 - 1/gamma**2)


def BetatoP(beta,rest_mass=938):
    """ returns p in [MeV/c]"""
    gamma = BetatoGamma(beta)
    return beta*gamma*rest_mass


def EtoP(E,rest_mass=938):
    """energy to impulsion"""
    p = sqrt(E**2 + 2*rest_mass*E)
    
    return p
  
def PtoE(p,rest_mass=938):
    """energy to impulsion"""
    E = sqrt(p**2 + rest_mass**2) - rest_mass
    
    return E          

def EtoR(E):
    R=exp(-0.013296*log(E)**3 + 0.15248*log(E)**2 + 1.2193*log(E) - 5.5064) 
    return R
  
def RtoE(R):
    E=exp(0.0015739*log(R)**3 - 0.0040274*log(R)**2 + 0.55919*log(R) + 3.4658)
    return E    

def EtoBrho(E,rest_mass=938):
    Brho  = 1/300*sqrt(E**2+2*rest_mass*E)
    return Brho

def PtoBrho(p):
    Brho  = 1/300*p
    return Brho

def Brho_scaling(refE,newE,rest_mass=938):
    refBrho = 1/rest_mass*sqrt(refE**2+2*rest_mass*refE)
    newBrho = 1/rest_mass*sqrt(newE**2+2*rest_mass*newE)
    
    return newBrho/refBrho
    

def Gaussian_fit(data,p0):
    """
    p0 is the initial guess for the fitting coefficients (A, mu and sigma above)# p0 is the initial guess for the fitting coefficients (A, mu and sigma below)
    """
    
    # Define model function to be used to fit to the data above:
    def gauss(x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    
    mask = ~np.isnan(data)
    
    hist, bin_edges = np.histogram(data[mask], density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    
    
    
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    
    print('Fitted amplitude = ', coeff[0])
    print('Fitted mean = ', coeff[1])
    print('Fitted standard deviation = ', coeff[2])
    
    return coeff




def EOM_particle_in_Bfield(Delta_z,beam,p_ref,Bx,By,Bz=0,part_mass=1.67*10**-27,rest_mass=938):
    #beam = np.array([z,sizeX,divX,sizeY,divY,0,dponp])
    
    
    
    
    light_speed = 2.998*10**8
    
    p_tot0 = p_ref * (beam[6] + 1)
    beta0 = PtoBeta(p_tot0,rest_mass)
    
    betaZ0 = beta0 / sqrt(1 + beam[2]**2 + beam[4]**2)
    betaX0 = betaZ0 * beam[2] 
    betaY0 = betaZ0 * beam[4]
    
    px0 = BetatoP(betaX0,rest_mass)
    py0 = BetatoP(betaY0,rest_mass)
    pz0 = BetatoP(betaZ0,rest_mass)
    
    
    
    Delta_t = Delta_z/betaZ0/light_speed
    
#    betaX0 = PtoBeta(px0,rest_mass)
#    betaY0 = PtoBeta(py0,rest_mass)
#    betaZ0 = PtoBeta(pz0,rest_mass)
    
    fact = 10**-6 * Delta_t * light_speed**2
    px = px0 + fact * (betaY0*Bz - betaZ0*By) 
    py = py0 + fact * (betaZ0*Bx - betaX0*Bz) 
    pz = pz0 + fact * (betaX0*By - betaY0*Bx) 
    
    # normalize to prevent propagation of rounding errors from linearization
    p_tot = sqrt(px**2 + py**2 + pz**2)
    px = px * p_tot0/p_tot
    py = py * p_tot0/p_tot
    pz = pz * p_tot0/p_tot
    
    
    
    
    # compute speeds
    vx = PtoV(px,rest_mass)
    vy = PtoV(py,rest_mass)
    vz = PtoV(pz,rest_mass)
    vx0 = betaX0*light_speed
    vy0 = betaY0*light_speed
    
    
    
    x = beam[1] + (vx+vx0)/2 * Delta_t 
    y = beam[3] + (vy+vy0)/2 * Delta_t 
    z = beam[0] + Delta_z
    
    
#    print('Delta_z = ',Delta_z)
#    print('Delta_t*vz = ',Delta_t*vz)
    
    divX = vx/vz
    divY = vy/vz
    
#    divX = (x-beam[1])/(z-beam[0])
#    divY = (y-beam[3])/(z-beam[0])
    
#    print('px = ',px)
#    print('py = ',py)
    
    
    dL = beam[5]
    dponp = beam[6]
    
    return [z,x,divX,y,divY,dL,dponp]



def gaussian(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:15:26 2019

@author: rvroerm
"""

import math as math
import numpy as np
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt


def transport_count_lines(transport_file,N_segments):
    # counts the number of relevant lines in file
    # N_segments is the number of points displayed by element
    
    it_z = 0    
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       cnt = 1
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           
           if data: # string is not empty
               if  data[0][0:2] == "4." or data[0][0:2] == "5." or data[0][0:3] == "50." or data[0][0:3] == "18.":
                   it_z = it_z + N_segments
                   print(data)
               elif data[0][0:2] == "3." or data[0][0:2] == "2.":
                   it_z = it_z + 1
                   print(data)
                   
           line = fp.readline()
           cnt += 1
           
    print('----')
    return it_z+1


def GTR_layout_from_transport(transport_file,coord,refE):
    
    
    
    angle = 0
    rot_mat = [[1,0],[0,1]]
    new_vect = [0,0]
    
    
    coord[0,:] = [0,0]
    it_z = 1
    
    Brho = 1/300*sqrt(refE**2+2*938*refE)
    
    filepath = transport_file
    with open(filepath) as fp:  
       line = fp.readline()
       cnt = 1
       
       while line:
           #print("Line {}: {}".format(cnt, line.strip()))
           
           data = line.split() 
           
           if data: # string is not empty
               
               if data[0][0:2] == "3.":
                   # drift
                   L = float(data[1].replace(";","") )
                   
                   coord[it_z,:] = coord[it_z-1,:] + np.matmul(rot_mat,[L,0])
                   
                   it_z = it_z + 1
                   
               if data[0][0:2] == "4.":
                   # sector bending
                   
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )
                   rho = Brho/B
                   B_angle = L/rho
                   
                   coord[it_z,:] = coord[it_z-1,:] + np.matmul(rot_mat,[rho*sin(B_angle),rho*(1-cos(B_angle))])
                   
                   angle = angle + B_angle
                   print('--')
                   print(angle)
                   print(B_angle)
                   
                   rot_mat = [[cos(angle),-sin(angle)],[sin(angle),cos(angle)]]
                   it_z = it_z + 1
                   
               if data[0][0:2] == "5.": 
                   # quad
                   L = float(data[1].replace(";","") )
                   new_vect[0] = L
                   new_vect[1] = 0
                   coord[it_z,:] = coord[it_z-1,:] + np.matmul(rot_mat,new_vect)
                   it_z = it_z + 1
                   
                   
           line = fp.readline()
           cnt += 1
       


    return coord  





def transport_input(transport_file,beam,refE,it_p,N_segments,gap,k1,k2,z,it_z):
    # opens a transport file and computes beam through magnets for one particle
    # The first particle it_p=0 should be the reference particle
    # all other particles are shown at a relative of that particle
    
    
    
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
                           [beam[:,:,it_p],it_z] = pole_face(angle,B,gap,k1,k2,beam[:,:,it_p],refE,it_z)
                           
                           
               if data[0][0:2] == "3.":
                   # drift
                   L = float(data[1].replace(";","") )
                   z = z + L
                   
                   [beam[:,:,it_p],it_z] = drift(L,beam[:,:,it_p],refE,it_z,1)
                   
               if data[0][0:2] == "4.":
                   # sector bending
                   
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )
                   n = float(data[3].replace(";","") )
                   
                   
                   if data_prev :
                       if data_prev[0] == "2.": # entrance pole face
                           angle = float(data_prev[1].replace(";","") )
                           [beam[:,:,it_p],it_z] = pole_face(angle,B,gap,k1,k2,beam[:,:,it_p],refE,it_z)
                           
                   
                   
                   [beam[:,:,it_p],it_z] = sec_bend(L,B,n,beam[:,:,it_p],refE,it_z,N_segments)
                   
                   
               if data[0][0:2] == "5.": 
                   # quad
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )
                   a = float(data[3].replace(";","") )/1000 # aperture converted to meters
                   
                   [beam[:,:,it_p],it_z] = quad(L,B,a,beam[:,:,it_p],it_z,N_segments,refE)
                   
               
               if data[0][0:3] == "18.": 
                   # sextupole
                   L = float(data[1].replace(";","") )
                   B = float(data[2].replace(";","") )
                   a = float(data[3].replace(";","") )/1000 # aperture converted to meters
                   angle = float(data[4].replace(";","") ) # not used in original transport
                   
                   [beam[:,:,it_p],it_z] = sextupole(L,B,a,angle,beam[:,:,it_p],it_z,N_segments,refE) 
                
                
               if data[0][0:3] == '50.':
                   # linear wedge
                   L = float(data[1].replace(";","") )
                   width = float(data[2].replace(";","") )
                   mat = data[3].replace(";","") 
                   
                   [beam[:,:,it_p],it_z] = wedge_X(L,width,mat,beam[:,:,it_p],it_z,N_segments,refE)
                
           prev_line = line     
           line = fp.readline()
           cnt += 1
    
    
    return [beam, it_z]




def drift(L,beam,refE,it_z,N_segments):
    
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



def RF_cavity(delta_E,cav_length,freq,phi,beam,refE,it_z,N_segments):
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
    

def quad(L,B,a,beam,it_z,N_segments,refE):
    # L =length, B=field, a=aperture
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    E = PtoE(p)
    gamma = PtoGamma(p,938)
    
    Brho  = 1/300*sqrt(E**2+2*938*E)
    #Brho  = 1/300*sqrt(refE**2+2*938*refE)
    k_q  = sqrt(abs(B)/a*(1/Brho))
    
    
    if B>0: 
        # X focusing
        quad_mat = np.array([[cos(k_q*L_segment),1/k_q*sin(k_q*L_segment),0,0,0,0],
                              [-k_q*sin(k_q*L_segment),cos(k_q*L_segment),0,0,0,0],
                            [0,0,cosh(k_q*L_segment),1/k_q*sinh(k_q*L_segment),0,0],
                            [0,0,k_q*sinh(k_q*L_segment),cosh(k_q*L_segment),0,0],
                            [0,0,0,0,1,L_segment/gamma**2],
                            [0,0,0,0,0,1]])
    else:
        # Y focusing
        quad_mat = np.array([[cosh(k_q*L_segment),1/k_q*sinh(k_q*L_segment),0,0,0,0],
                              [k_q*sinh(k_q*L_segment),cosh(k_q*L_segment),0,0,0,0],
                            [0,0,cos(k_q*L_segment),1/k_q*sin(k_q*L_segment),0,0],
                            [0,0,-k_q*sin(k_q*L_segment),cos(k_q*L_segment),0,0],
                            [0,0,0,0,1,L_segment/gamma**2],
                            [0,0,0,0,0,1]])
    
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + L_segment
        beam[it_z,1:7] = np.matmul(quad_mat,np.transpose(beam[it_z-1,1:7]))
        
        L_eff = abs(L_segment/cos(sqrt(beam[it_z,2]**2+beam[it_z,4]**2))) # effective length travelled by particle (negliect change in direction)
        beam[it_z,5] = (L_segment - L_eff) + beam[it_z,5] #higher order correction related to beam angle
    
    
    return [beam, it_z]


def sextupole(L,B,a,alpha,beam,it_z,N_segments,refE):
    # sextupole: L =length, B=field, a=aperture, alpha = angle
    
    alpha = math.radians(alpha)
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    E = PtoE(p)
    gamma = PtoGamma(p,938)
    
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

def sec_bend(L,B,n,beam,refE,it_z,N_segments):
    # NEED TO CLARIFY DIFFERENCE IN EFFECTIVE LEGNTH
    
    
#    L_segment = L/N_segments
#    
#    Brho  = 1/300*sqrt(refE**2+2*938*refE)
    
    
    
    L_segment = L/N_segments
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    E = PtoE(p)
    gamma = PtoGamma(p,938)
    
    Brho  = 1/300*sqrt(E**2+2*938*E)
    
    
    h = B/Brho
    kx  = sqrt((1-n)*h**2)
    ky  = sqrt(n*h**2) 
    
    
    if ky!=0:
        bend_mat = np.array([[cos(kx*L_segment),1/kx*sin(kx*L_segment),0,0,0,h/kx**2*(1-cos(kx*L_segment))],
                              [-kx*sin(kx*L_segment),cos(kx*L_segment),0,0,0,h/kx*sin(kx*L_segment)],
                              [0,0,cos(ky*L_segment),1/ky*sin(ky*L_segment),0,0],
                              [0,0,-ky*sin(ky*L_segment),cos(ky*L_segment),0,0],
                              [-h/kx*sin(kx*L_segment),-h/kx**2*(1-cos(kx*L_segment)),0,0,1,-h**2/kx**3*(kx*L_segment-sin(kx*L_segment))],
                              [0,0,0,0,0,1]])
    else:
        bend_mat = np.array([[cos(kx*L_segment),1/kx*sin(kx*L_segment),0,0,0,h/kx**2*(1-cos(kx*L_segment))],
                              [-kx*sin(kx*L_segment),cos(kx*L_segment),0,0,0,h/kx*sin(kx*L_segment)],
                              [0,0,1,L_segment,0,0],
                              [0,0,0,1,0,0],
                              [-h/kx*sin(kx*L_segment),-h/kx**2*(1-cos(kx*L_segment)),0,0,1,-h**2/kx**3*(kx*L_segment-sin(kx*L_segment))],
                              [0,0,0,0,0,1]])
    
    
    for i in range(0,N_segments):
        it_z = it_z + 1
        beam[it_z,0] = beam[it_z-1,0] + L_segment
        beam[it_z,1:7] = np.matmul(bend_mat,np.transpose(beam[it_z-1,1:7]))
    
    
    return [beam, it_z]


def pole_face(angle,B,gap,k1,k2,beam,refE,it_z):
    # k1 and k2 related to frindge field properties; default k1=0.5, k2=0
    # see SLAC transport manual
    
    
    ref_p = EtoP(refE)
    p = ref_p + beam[it_z,6]*ref_p
    E = PtoE(p)
    
    Brho  = 1/300*sqrt(E**2+2*938*E)
    
    radius = Brho/B
    
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




def bending(L,B,n,pole_in,pole_out,gap,k1,k2,beam,refE,it_z,N_segments):
    
    [beam, it_z] = pole_face(pole_in,B,gap,k1,k2,beam,refE,it_z)
    
    [beam, it_z] = sec_bend(L,B,n,beam,refE,it_z,N_segments)
    
    [beam, it_z]= pole_face(pole_out,B,gap,k1,k2,beam,refE,it_z)

    return [beam, it_z]





def solenoid(L,B,a,beam,it_z,N_segments,refE):
    
    
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


    


def wedge_X(thickness_max,width,mat,beam,it_z,N_segments,refE):
    
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


def collimator(L_m,a,mat,beam,it_z,N_segments,refE):
    # L_m=collimator thickness, a=aperture radius, mat=material
    
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
    # computes scattering angle

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
            
            beam_out[5] = (new_p-ref_p)/ref_p
            
        else:
            beam_out = np.empty(6)
            beam_out[:] = np.nan
            sigma = 0
             
    

    return [beam_out, sigma]


def PtoGamma(p,rest_mass):
    # impulsion (in Mev/c) to gamma factor
    
    beta = p/sqrt(rest_mass**2 + p**2)
    gamma = 1/sqrt(1-beta**2)
    
    return gamma


def PtoV(p,rest_mass):
    # impulsion (in Mev/c) to speed (in m/s)
    
    light_speed = 3*10**8
    speed = p/sqrt(rest_mass**2 + p**2)*light_speed
    
    return speed
    

def EtoP(E):
    #energy to impulsion
    p = sqrt(E**2 + 2*938*E)
    
    return p
  
def PtoE(p):
    #energy to impulsion
    E = sqrt(p**2 + 938**2) - 938
    
    return E          

def EtoR(E):
    R=exp(-0.013296*log(E)**3 + 0.15248*log(E)**2 + 1.2193*log(E) - 5.5064) 
    return R
  
def RtoE(R):
    E=exp(0.0015739*log(R)**3 - 0.0040274*log(R)**2 + 0.55919*log(R) + 3.4658)
    return E    

def EtoBrho(E):
    Brho  = 1/300*sqrt(E**2+2*938*E)
    return Brho


def Gaussian_fit(data,p0):
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)# p0 is the initial guess for the fitting coefficients (A, mu and sigma below)
    
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



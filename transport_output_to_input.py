# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:40:18 2019

@author: rvroerm
"""

import math as math
import numpy as np
import os
from math import sin, cos, tan, sinh, cosh, tanh, exp, log, log10, sqrt
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt




input_file = "C:/TRANS/for001.dat"
output_file = "D:/temp/transport_output.txt" #copy Transport output in that file

temp_file = "C:/TRANS/temp.txt"


with open(input_file) as fp1,  open(output_file) as fp2 , open(temp_file,'w') as fp3: 
    line_input = fp1.readline()
    cnt_input = 1
    
    line_output = fp2.readline()
    cnt_output = 1
    
   
    while line_input:
       
       data_input = line_input.split() 
       
       if data_input: # string is not empty
           if data_input[0][0:2] == '5.':
               
               while line_output:
                   data_output = line_output.split() 
                   if data_output:
                       if data_output[0] == '*QUAD*':
                           line_output = fp2.readline()
                           cnt_output += 1
                           break # found corresponding quad in output file
                   line_output = fp2.readline()
                   cnt_output += 1
                   
               if data_output[6] == 'T':
                   # no label
                   data_input[1] = str(float(data_output[1]))
                   data_input[3] = str(float(data_output[3]))
                   data_input[2] = str(float(data_output[5]))
               else:
                   data_input[1] = str(float(data_output[2]))
                   data_input[3] = str(float(data_output[4]))
                   data_input[2] = str(float(data_output[6]))
                   
               fp3.write(' '.join(data_input) + '\n')
               
           elif data_input[0][0:2] == '4.':
               while line_output:
                   # Replace the target string (correct bug from Transport)
                   line_output = line_output.replace('T', 'T ')
                   
                   data_output = line_output.split() 
                   if data_output:
                       if data_output[0] == '*BEND*':
                           line_output = fp2.readline()
                           cnt_output += 1
                           break # found corresponding quad in output file
                   line_output = fp2.readline()
                   cnt_output += 1
               
               if data_output[9] == 'Deg':
                   # no label
                   data_input[1] = str(float(data_output[1]))
                   data_input[2] = str(float(data_output[5]))
                   data_input[3] = str(float(data_output[7]))
               else:
                   data_input[1] = str(float(data_output[2]))
                   data_input[2] = str(float(data_output[6]))
                   data_input[3] = str(float(data_output[8])) 
                   
               fp3.write(' '.join(data_input) + '\n')
               
           else:
               fp3.write(' '.join(data_input) + '\n')
       else:
           fp3.write('\n')  
       
        
       line_input = fp1.readline()    
       cnt_input += 1

       
fp2.close()
fp1.close() 
fp3.close()


# replace initial file with temp file

os.remove(input_file)
os.rename(temp_file, input_file)
       
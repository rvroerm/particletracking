# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 12:18:28 2019

@author: rvroerm
"""



from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt

def funcLine(param,x):
    return param[0]*x+param[1]

def errorFuncLine(param,x,y):
    return funcLine(param,x)-y

def funcQuad(param,x):
    return param[0]*x**2+param[1]*x+param[2]

def errorfuncQuad(param,x,y):
    return funcQuad(param,x)-y


def funcGaussian(param,x):
    return param[0] * np.exp( - (x - param[1])**2.0 / (2.0 * param[2]**2.0) ) 

def errorFuncGaussian(param,x,y):
    return funcGaussian(param,x)-y

def funcDoubleGaussian(param,x):
    return param[0] * np.exp( - (x - param[1])**2.0 / (2.0 * param[2]**2.0) ) \
         + param[3] * np.exp( - (x - param[4])**2.0 / (2.0 * param[5]**2.0) ) 

def errorFuncDoubleGaussian(param,x,y):
    return funcDoubleGaussian(param,x)-y

def main():
   # data provided
   x=np.array([1.0,2.5,3.5,4.0,1.1,1.8,2.2,3.7])
   y=np.array([6.008,15.722,27.130,33.772,5.257,9.549,11.098,28.828])
   
   
   # func is going to be a placeholder for funcLine,funcQuad or whatever 
   # function we would like to fit
   func=funcGaussian
   
   # ErrorFunc is the diference between the func and the y "experimental" data
   ErrorFunc=lambda param,x,y: func(param,x)-y
   
   
   #tplInitial contains the "first guess" of the parameters 
   param_initial=(100,0,0.003)
   # leastsq finds the set of parameters in the tuple tpl that minimizes
   # ErrorFunc=yfit-yExperimental
   param_final,success=leastsq(ErrorFunc,param_initial[:],args=(x,y))
   print(" gaussian fit ",param_final)
   xx1=np.linspace(x.min(),x.max(),50)
   yy1=func(param_final,xx1)
   
   plt.figure(8)
   plt.plot(x,y,'b')
   plt.plot(xx1,yy1,'r-')
   plt.show()

if __name__=="__main__":
   main()

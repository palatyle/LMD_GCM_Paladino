#! /usr/bin/env python
#-*- coding:Utf-8 -*-

from    netCDF4               import    Dataset
from	scipy.ndimage.filters                import	gaussian_filter1d	
from	scipy.ndimage.filters                import	convolve1d
from	scipy.ndimage.filters                import	uniform_filter1d		
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
from math import log

#********** INPUT PARAMETERS **********

data1=np.loadtxt('bands_sun.txt')



        
mpl.plot((10000./data1[:,0]+10000./data1[:,1])/2.,data1[:,2]/((10000./data1[:,0]-10000./data1[:,1])),'bs')        
mpl.plot((10000./data1[:,0]+10000./data1[:,1])/2.,data1[:,2]/((10000./data1[:,0]-10000./data1[:,1])))           
        
        
for i in range(0,36,1) :       
   print (10000./data1[i,0]+10000./data1[i,1])/2., data1[i,2]/((10000./data1[i,0]-10000./data1[i,1]))
        
        
       

mpl.show()





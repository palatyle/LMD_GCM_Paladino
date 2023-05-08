#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from numpy import *                
from    netCDF4               import    Dataset
from	scipy.ndimage.filters                import	gaussian_filter1d	
from	scipy.ndimage.filters                import	convolve1d
from	scipy.ndimage.filters                import	uniform_filter1d		
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
from math import log
from matplotlib.colors import LogNorm

#********** INPUT PARAMETERS **********

g = loadtxt('g.dat')

data1 = loadtxt('corrk_gcm_IR_mod.dat')

wn = loadtxt('narrowbands_IR.in')

band=38 # Number of bands
ngauss=17 # Number of Gauss points
temp=10 # Number of temperature points
pres=9 # Number of pressure points
x_vmr=11 # Number of mixing ratio points

corrk=np.zeros((temp,pres,x_vmr,band,ngauss),dtype='f')

count=0
for i in range(0,ngauss,1):
   for j in range(0,band,1):
      for k in range(0,x_vmr,1):
         for l in range(0,pres,1):
            for m in range(0,temp,1):
               corrk[m,l,k,j,i]=data1[count]
               count=count+1

print g

g_cum=np.zeros(len(g),dtype='f')

g_cum[0]=g[0]
for i in range(1,len(g),1):
   g_cum[i]=g_cum[i-1]+g[i]

print g_cum

wn_mean=np.zeros(len(wn),dtype='f')

for i in range(0,len(wn),1):
   wn_mean[i]=(wn[i,0]+wn[i,1])/2.
   
print wn_mean

mpl.figure(1)
mpl.ylabel('absorption', fontsize=15)
mpl.xlabel('distribution', fontsize=15)
mpl.semilogy(g_cum,corrk[0,0,0,0,:])

mpl.show()

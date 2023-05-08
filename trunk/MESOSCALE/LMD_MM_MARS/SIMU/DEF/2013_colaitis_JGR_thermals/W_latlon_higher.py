#! /usr/bin/env python

from netCDF4 import Dataset
import matplotlib.pyplot as mpl
import numpy as np
from os import system
from matplotlib.pyplot import contour,contourf,colorbar
from matplotlib.cm import get_cmap
from myplot import makeplotres

folder='/san/work/colaitis/SIMUS/LES.833x833x133.CaseC_zipbl_64proc/LMD_LES_MARS_LARGE_DOMAIN.23473/'
file='wrfout_d01_9999-01-01_09:15:00'
path=folder+file
nc=Dataset(path)
varW=nc.variables["W"][0,60,:,:]
dimensions=np.array(varW).shape
coordx=np.arange(dimensions[0])
coordy=np.arange(dimensions[1])

min=-8.
max=8.
zelevels = np.linspace(min,max,num=30)
contourf(coordx, coordy, varW, zelevels, cmap = get_cmap(name='hsv'))
colorbar(fraction=0.05,pad=0.03,format='%.1f',orientation='vertical',\
                                      ticks=np.linspace(min,max,num=30),extend='neither',spacing='proportional' )
mpl.title('Vertical velocity at level 60')
mpl.xlabel('East-West')
mpl.ylabel('South-North')
makeplotres('W_latlon_60',res=200.,disp=False)

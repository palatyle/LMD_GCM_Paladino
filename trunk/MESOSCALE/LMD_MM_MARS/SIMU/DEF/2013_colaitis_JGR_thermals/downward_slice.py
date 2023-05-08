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
varQ=nc.variables["downward"][0,:,0,:]
dimensions=np.array(varQ).shape
print 'dimensions',dimensions
coordx=np.arange(dimensions[0])
coordy=np.arange(dimensions[1])

min=0.
max=0.00005
zelevels = np.linspace(min,max,num=30)
contourf(coordy, coordx, varQ, zelevels, cmap = get_cmap(name='jet'))
colorbar(fraction=0.05,pad=0.1,format='%.1e',orientation='horizontal',\
                                      ticks=np.linspace(min,max,num=6),extend='neither',spacing='proportional' )
##mpl.axis('equal')
mpl.xlim(xmin=0)
mpl.xlim(xmax=832)
mpl.ylim(ymin=0)
mpl.ylim(ymax=132)
mpl.xlabel('East-West')
mpl.ylabel('Bottom-Top')
mpl.axes().set_aspect('equal')
mpl.title('downward tracer concentration, East-West slice')
makeplotres('downward_slice',res=200.,disp=False)

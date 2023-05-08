#! /usr/bin/env python

import numpy as np
from math import *
from struct import *
from netCDF4 import Dataset
from array import *

filename = 'Relief.nc'
file=Dataset(filename,mode='r')
topo=file.variables['RELIEF'][:]
res=len(topo)
topo=np.flipud(topo)
topo=topo*1. + 9000.

part2b=[]
partb=[]
for j in range(res) :
  a=[]
  b=[]
  for i in range(res) :
     a.append(pack('>h',topo[j][2*res-1-i]))
     b.append(pack('>h',topo[j][res-i]))
  partb.append(a)
  part2b.append(b) 

       # Eastern part
#f = open('00181-00360.00001-00180', 'wb') # for low res : 180x360 pts
#f = open('04097-08192.00001-04096', 'wb') # for high res : 4096x8192 pts
f.write(np.array(part2b))
	# Western part
#f = open('00001-00180.00001-00180','wb')
#f = open('00001-04096.00001-04096','wb')
f.write(np.array(partb))

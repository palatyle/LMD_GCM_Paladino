#! /usr/bin/env python

import numpy as np
from math import *
from struct import *
from netCDF4 import Dataset
from array import *

#filename = 'HR_Relief.nc'
#file=Dataset(filename,mode='r')
#topo=file.variables['RELIEF'][:]

#res=len(topo)
res=180
#topo=topo*1. + 9000.
topo=0
topo=topo + 9000.

part2b=[]
partb=[]
for j in range(res) :
  a=[]
  b=[]
  for i in range(res) :
     a.append(pack('>h',topo))
     #a.append(pack('>h',topo[j][i]))
     #a.append(pack('>h',topo[j][2*res-1-i]))
     b.append(pack('>h',topo))
     #b.append(pack('>h',topo[j][i+res]))
     #b.append(pack('>h',topo[j][res-i]))
  partb.append(a)
  part2b.append(b)
       # Eastern part
#f = open('00res-00res*2.00001-00res', 'wb')
f = open('00181-00360.00001-00180', 'wb')
f.write(np.array(part2b))
	# Western part
#f = open('00001-00res.00001-00res','wb')
f = open('00001-00180.00001-00180','wb')
f.write(np.array(partb))

#! /usr/bin/env python

##file:///home/aymeric/Software/epd-7.0-2-rh5-x86/Doc/library/stdtypes.html?highlight=file#file.write
##file:///home/aymeric/Software/epd-7.0-2-rh5-x86/Doc/library/struct.html?highlight=endian

## marche pas probleme encodage...

from myplot import getfield,definesubplot,smooth
from netCDF4 import Dataset
from matplotlib.pyplot import contourf,show,pcolor,subplot,figure
from numpy import flipud,array,transpose
from os import sys

print sys.byteorder

charvar = "albedo"
charvar = "zMOL"
#charvar = "thermal"
charvar = "z0"

nc = Dataset("surface.nc")

var = getfield(nc,charvar)*1000.
var = flipud(var)
var = smooth(var,10) ##change le type!

tile = 180  # resolution

print var.shape

# Eastern part
epart = var[0:tile-1,tile:2*tile-1]
epart = transpose(epart)
fid = open('00181-00360.00001-00180','wb')
#fid.write(str(part2))
#pickle.dump(epart,fid)
#epart = array(epart,'<h')
epart = array(epart,dtype='h')
epart.tofile(fid)
fid.close()
#### integer*2
print epart.itemsize

## Western part
wpart = var[0:tile-1,0:tile-1]
wpart = transpose(wpart)
fid2 = open('00001-00180.00001-00180','wb')
#fid2.write(str(part))
#pickle.dump(wpart,fid2)
#wpart = array(wpart,'<h')
wpart = array(wpart,dtype='h')
wpart.tofile(fid2)
fid2.close()

#fid = open('00181-00360.00001-00180','r')
#epart = pickle.load(fid)
#fid.close()
#
#fid2 = open('00001-00180.00001-00180','r')
#wpart = pickle.load(fid2)
#fid2.close()

fig = figure()
sub = definesubplot(2,fig)

subplot(sub)
contourf(epart)

subplot(sub+1)
contourf(wpart)

show()



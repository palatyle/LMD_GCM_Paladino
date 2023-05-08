#!/usr/bin/env python
# coding: utf-8

import argparse
import netCDF4

## parse script arguments
parser=argparse.ArgumentParser(description="Load a variable from a first file to (over-)write it in a second file")
## add a mandatory argument
parser.add_argument("file1",help="first input file name")
## add another mandatory argument
parser.add_argument("file2",help="second input/output file name")
## add another mandatory argument
parser.add_argument("variable",help="variable to read from first file and (over-)write to second file")

# get arguments
args=parser.parse_args()
filename1=args.file1
filename2=args.file2
varname=args.variable

print "Replacing variable ",varname," in ",filename2," by the one in ",filename1

# open files
try:
  dset1=netCDF4.Dataset(filename1,"r") # read only
except:
  print "Error: cannot find file ",filename1
  exit() 

try:
  dset2=netCDF4.Dataset(filename2,"a") # read/write
except:
  print "Error: cannot find file ",filename2
  exit() 

# load data
try:
  data1=dset1.variables[varname]
except:
  print "Error, cannot find variable ",varname," in file ",filename1
  exit() 

try:
  data2=dset2.variables[varname]
except:
  print "Error, cannot find variable ",varname," in file ",filename2
  exit() 

# copy data1 into data2
if data2.shape != data1.shape:
  print "Error: shape mismatch for variable ",varname," between the 2 files!"
  print "    ",data1.shape," in ",filename1
  print "    ",data2.shape," in ",filename2
  exit()
nd2 = data2.ndim # number of dimensions of data2
if   nd2 == 1: data2[:]=data1[:]
elif nd2 == 2: data2[:,:]=data1[:,:]
elif nd2 == 3: data2[:,:,:]=data1[:,:,:]
elif nd2 == 4: data2[:,:,:,:]=data1[:,:,:,:]

# close file to write data to file
dset2.close()

print "All's well that ends well."


#! /usr/bin/octave -qf

# -----------------------------------------------
# thermal_build.m					
#   Script to be used with GNU-Octave
#   > For NETCDF support, get 'octave-forge' pckg octcdf (install headers)
#   and install it in octave with 'pkg install octcdf-1.0.11.tar.gz'
#   (you also need libnetcdf-dev)
# -----------------------------------------------
# Purpose:
#   Thermal inertia GCM NETCDF input
#	>> WRF geogrid tiles
# Author:
#   A. Spiga - 03/2007
# -----------------------------------------------

# Locate NETCDF file
filename = 'surface.nc';
resolution = 180;

# Open NETCDF file - read-only
nc = netcdf('surface.nc','r')

# Read ...
# Thermal inertia
nv = ncvar(nc){4};
thermal = nv(:);
size(thermal)

# Flip North/South 
thermal = flipud(thermal);
#thermal = fliplr(thermal);

# Create 2 WRF data tiles for geogrid
tile=resolution;
	# Eastern part
	part = thermal(1:1:tile,1:1:tile)';
	#fid = fopen('00181-00360.00001-00180','wb','b');
	fid = fopen('00001-00180.00001-00180','wb','b');
	fwrite(fid,part,'integer*2');

	# Western part
	part2 = thermal(1:1:tile,tile+1:1:2*tile)';
	#fid = fopen('00001-00180.00001-00180','wb','b');
	fid = fopen('00181-00360.00001-00180','wb','b');
	fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:1:tile,1:1:tile);
yeah2 = part2(1:1:tile,1:1:tile);
contour(yeah)
contour(yeah2)

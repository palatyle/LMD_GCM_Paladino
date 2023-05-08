#! /usr/bin/octave -qf

# -----------------------------------------------
# tesalbedo_build.m					
#   Script to be used with Matlab or GNU-Octave
# -----------------------------------------------
# Purpose:
#   MOLA MEGDR binary file >> WRF geogrid tiles
# Author:
#   A. Spiga - 03/2007
# -----------------------------------------------

# Locate MOLA binary file
filename = 'global_albedo_8ppd.img';
resolution = 8;

# Read topographical data (PC_REAL, 32-bits/4-bytes float)
f = fopen(filename,'r','ieee-le');
el = fread(f,[360*resolution Inf],'float32')';

# Flip North/South 
el = flipud(el);
# Scale factor (ie accuracy): 10000
el = 10000.*el;
# Conversion float >> integer
el = round(el);

# Create 2 WRF data tiles for geogrid
tile=180*resolution;
	# Eastern part
	part = el(1:1:tile,1:1:tile)';
	fid = fopen('00001-01440.00001-01440','wb','b');
	fwrite(fid,part,'integer*2');

	# Western part
	part2 = el(1:1:tile,tile+1:1:2*tile)';
	fid = fopen('01441-02880.00001-01440','wb','b');
	fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:10:tile,1:10:tile);
yeah2 = part2(1:10:tile,1:10:tile);
contour(yeah)
contour(yeah2)


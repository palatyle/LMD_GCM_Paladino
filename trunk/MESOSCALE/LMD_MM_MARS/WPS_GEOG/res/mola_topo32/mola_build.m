#! /usr/bin/octave -qf


# -----------------------------------------------
# mola_build.m					
#   Script to be used with Matlab or GNU-Octave
# -----------------------------------------------
# Purpose:
#   MOLA MEGDR binary file >> WRF geogrid tiles
# Author:
#   A. Spiga - 03/2007
# -----------------------------------------------

# Locate MOLA binary file
filename = 'megt90n000fb.img';
resolution = 32;

# Read topographical data (MSB/big endian, 16-bits/2-bytes integer)
f = fopen(filename,'r','ieee-be');
el = fread(f,[360*resolution Inf],'int16')';

# Get rid of negative values and flip North/South 
el = el + 9000;
el = flipud(el);

# Create 2 WRF data tiles for geogrid
tile=180*resolution;
	# Eastern part
	part = el(1:1:tile,1:1:tile)';
	fid = fopen('05761-11520.00001-05760','wb','b');
	fwrite(fid,part,'integer*2');

	# Western part
	part2 = el(1:1:tile,tile+1:1:2*tile)';
	fid = fopen('00001-05760.00001-05760','wb','b');
	fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:100:tile,1:100:tile);
yeah2 = part2(1:100:tile,1:100:tile);
contour(yeah)
contour(yeah2)

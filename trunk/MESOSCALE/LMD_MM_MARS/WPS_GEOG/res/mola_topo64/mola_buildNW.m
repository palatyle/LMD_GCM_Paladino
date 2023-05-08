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
filename = 'megt90n000gb.img';
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
#
#just have to change the names !
#-- xrange - yrange
	fid = fopen('11521-17280.05761-11520','wb','b');
	fwrite(fid,part,'integer*2');

	# Western part
	part2 = el(1:1:tile,tile+1:1:2*tile)';
	fid = fopen('17281-23040.05761-11520','wb','b');
	fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:300:tile,1:300:tile);
yeah2 = part2(1:300:tile,1:300:tile);
contour(yeah)
contour(yeah2)

#! /usr/bin/octave -qf

# -----------------------------------------------------
# paleo_build.m
#   Script to be used with Matlab or GNU-Octave
# -----------------------------------------------------
# Purpose:
#   MOLA_like MEGDR binary file >> WRF geogrid tiles
# Author:
#   A. Spiga - 05/2011
# -----------------------------------------------------

# Locate MOLA binary file
filename = 'bu_mola64_v7.flt';
resolution = 64;

# Read topographical data 
###32-bit signed integer, ieee floating point format, big endian, 64 pixels per degree
###...in reality this is different...
f = fopen(filename,'r','ieee-le');
el = fread(f,[23040 5760],'float32')';

# Get rid of negative values and flip North/South
el = el + 9000;
el = flipud(el);

# Create 2 WRF data tiles for geogrid
tile=180*resolution;
tilered=180*resolution/2;

        # Eastern part
        part = el(1:1:tilered,1:1:tile)';
        #fid = fopen('11521-23040.00001-05760','wb','b');
        fid = fopen('11521-23040.05761-11520','wb','b');
        fwrite(fid,part,'integer*2');

        # Western part
        part2 = el(1:1:tilered,tile+1:1:2*tile)';
        #fid = fopen('00001-11520.00001-05760','wb','b');
        fid = fopen('00001-11520.05761-11520','wb','b');
        fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:100:tile,1:100:tilered);
yeah2 = part2(1:100:tile,1:100:tilered);
contour(yeah)
contour(yeah2)

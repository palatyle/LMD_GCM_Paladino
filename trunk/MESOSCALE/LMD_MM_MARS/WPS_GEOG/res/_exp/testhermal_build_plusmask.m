# -----------------------------------------------
# testhermal_build.m					
#   Script to be used with Matlab or GNU-Octave
# -----------------------------------------------
# Purpose:
#   MOLA MEGDR binary file >> WRF geogrid tiles
# Author:
#   A. Spiga - 03/2007
# -----------------------------------------------

# Locate MOLA binary file
filename = 'NBmap2007.bin';
filename2 = 'NBmsk2007.bin';
resolution = 20;

# Read topographical data (SUN, 16-bits/2-bytes integer)
f = fopen(filename,'r','ieee-be');
	bintitle = fread(f, 14400, 'char');   %read in the header
	title = char(bintitle');
el = fread(f,[360*resolution Inf],'int16')';

# Read mask data (SUN, 8-bits/1-byte integer)
f2 = fopen(filename2,'r','ieee-be');
	bintitle = fread(f2, 7200, 'char');   %read in the header
	title = char(bintitle');
el2 = fread(f2,[360*resolution Inf],'int8')';

# Apply mask by multiplying term-by-term the two arrays
# NB: all missing values will be 0.
el = el .* el2;


## Flip East/West 
el = fliplr(el);
## Scale factor (ie accuracy): 10000
#el = 10000.*el;
## Conversion float >> integer
#el = round(el);

# Create 2 WRF data tiles for geogrid
tile=180*resolution;
	# Eastern part
	part = el(1:1:tile,1:1:tile)';
	fid = fopen('03601-07200.00001-03600','wb','b');
	fwrite(fid,part,'integer*2');

	# Western part
	part2 = el(1:1:tile,tile+1:1:2*tile)';
	fid = fopen('00001-03600.00001-03600','wb','b');
	fwrite(fid,part2,'integer*2');

# Check the resulting arrays
yeah = part(1:100:tile,1:100:tile);
yeah2 = part2(1:100:tile,1:100:tile);
contour(yeah)
contour(yeah2)

#disp(yeah2)


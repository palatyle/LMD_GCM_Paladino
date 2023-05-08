# -----------------------------------------------
# albedo_build.m					
#   Script to be used with GNU-Octave
#   > For NETCDF support, get 'octave-forge' pckg 
# -----------------------------------------------
# Purpose:
#   Surface albedo GCM NETCDF input
#	>> WRF geogrid tiles
# Author:
#   A. Spiga - 03/2007
# -----------------------------------------------

# Locate NETCDF file
filename = 'albedo.nc';
resolution = 501;

# Open NETCDF file - read-only
nc = netcdf('albedo.nc','r')

# Read ...
# Surface albedo
nv = ncvar(nc){3};
albedo = nv(:);
size(albedo)

# Flip North/South 
albedo = flipud(albedo);
albedo = fliplr(albedo);

# Scale factor: 100
albedo = 100.*albedo;


# Create 2 WRF data tiles for geogrid
tile=resolution;

	# Western part
	part = albedo(1:1:tile,1:1:tile)';
	fid = fopen('00001-00501.00001-00501','wb','b');
	fwrite(fid,part,'integer*2');

# Check the resulting arrays
yeah = part(1:1:tile,1:1:tile);
contour(yeah)


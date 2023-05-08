% -----------------------------------------------
% custom_build.m					
%   Script to be used with Matlab or GNU-Octave
% -----------------------------------------------
% Purpose:
%   Custom topography in a peculiar region >> WRF geogrid tiles
% Author:
%   A. Spiga - 03/2007
% -----------------------------------------------

% Locate ENVI data file
filename = 'Large.dat';

% Read topographical data (PC linux, 16-bits/2-bytes integer)
f = fopen(filename,'r','ieee-le');
el = fread(f,[641 Inf],'int16')';

% Topographical offset to manage negative values 
el = el + 9000.;

% Create WRF data tiles (with 0. filling, to achieve a square tile)
tilex=641;
tiley=385;
part = zeros(tilex,tilex);
part(1:tilex,1:tiley) = el(1:1:tiley,1:1:tilex)';
fid = fopen('00001-00641.00001-00641','wb','b');
fwrite(fid,part,'integer*2');

% Check the resulting arrays
yeah = part(1:10:tilex,1:10:tilex);
contour(yeah)



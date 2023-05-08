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
filename = '4000.dat';

% Read topographical data (PC linux, 32-bits/4-bytes integer)
f = fopen(filename,'r','ieee-le');
el = fread(f,[1000 Inf],'float32')';

%% Scale factor (ie accuracy): 100
%el = 100.*el;
%% Topographical offset to manage negative values
%el = el + 9000.*100.;
% Conversion float >> integer
el = round(el);
el = el + 9000.;

% Create WRF data tiles (with 0. filling, to achieve a square tile)
tilex=1000;
%tiley=1000;
%part = zeros(tilex,tilex);
part(1:1:tilex,1:1:tilex) = el(1:1:tilex,1:1:tilex)';
fid = fopen('00001-01000.00001-01000','wb','b');
fwrite(fid,part,'integer*2');

% Check the resulting arrays
yeah = part(1:10:tilex,1:10:tilex);
contour(yeah)



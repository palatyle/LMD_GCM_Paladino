% -----------------------------------------------
% custom_build.m
%   Script to be used with Matlab or GNU-Octave
% -----------------------------------------------
% Purpose:
%   Custom topography in a peculiar region >> WRF geogrid tiles
% Author:
%   A. Spiga - 03/2007
% -----------------------------------------------

load 'tyler.mat';
tilex=320
tiley=5760
el = alb(1:1:tilex,1:1:tiley)';

% Scale factor (ie accuracy): 10000
el = 10000.*el;
% Conversion float >> integer
el = round(el);

disp("coordinates")
disp(LON(1,1:500:tiley))
disp(LAT(1:20:tilex,1)')
disp(LON(1,2)-LON(1,1))
disp(LAT(2,1)-LAT(1,1))


yeah = el(2000,1:10:tilex);
plot(yeah)

%% Create WRF data tiles (with 0. filling, to achieve a square tile)
%tilex=641;
%tiley=385;

tilexx=5760
tileyy=2880

part = zeros(tilexx,tileyy);
part(1:tilexx,tileyy-tilex+1:tileyy) = el;
fid = fopen('00001-05760.00001-02880','wb','b');
fwrite(fid,part,'integer*2');

plot(part(2000,1:tileyy))

plot(part(1000,1:tileyy))


%sho=part(1:100:tiley,1:100:tiley);
%contour(sho)


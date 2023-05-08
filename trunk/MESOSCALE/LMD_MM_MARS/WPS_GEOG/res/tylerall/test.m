% -----------------------------------------------
% custom_build.m
%   Script to be used with Matlab or GNU-Octave
% -----------------------------------------------
% Purpose:
%   Custom topography in a peculiar region >> WRF geogrid tiles
% Author:
%   A. Spiga - 03/2007
% -----------------------------------------------

load 'TylerData_8ppd.mat';
tilex=1440
tiley=2880
el = alb(1:1:tilex,1:1:tiley)';

% Scale factor (ie accuracy): 10000
el = 10000.*el;
% Conversion float >> integer
el = round(el);

disp("coordinates")
disp(lon(1,1:200:tiley))
disp(lat(1:200:tilex,1)')
disp(lon(1,2)-lon(1,1))
disp(lat(2,1)-lat(1,1))


yeah = el(2000,1:10:tilex);
plot(yeah)

%% Create WRF data tiles (with 0. filling, to achieve a square tile)
%tilex=641;
%tiley=385;

tilexx=2880
tileyy=1440

part = zeros(tilexx,tileyy);
part(1:tilexx,1:tileyy) = el;
fid = fopen('albedo/00001-02880.00001-01440','wb','b');
fwrite(fid,part,'integer*2');

plot(part(2000,1:tileyy))

plot(part(1000,1:tileyy))


%% now thermal inertia
el = ti(1:1:tilex,1:1:tiley)';
part = zeros(tilexx,tileyy);
part(1:tilexx,1:tileyy) = el;
fid = fopen('ti/00001-02880.00001-01440','wb','b');
fwrite(fid,part,'integer*2');


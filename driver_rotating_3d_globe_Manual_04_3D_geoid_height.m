% m-file: driver_rotating_3d_globe_Manual_04_3D_geoid_height
% Examples in the Manual: 4. Geoid height in 3D image/animation
% folder: private/fig04
%
% When using this driver or fuctions it calls, please cite the reference:
%  Bezdek A, Sebera J, 2013. MATLAB script for 3D visualizing geodata on a rotating globe.
%  Computers & geosciences 56, 127�130. http://dx.doi.org/10.1016/j.cageo.2013.03.007

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015
% 6/2016 update to new Matlab versions up to R2015a

%% Selection of geopotential model and computation/loading of the grid values
% model='asu-ch-0309'; %nmax 120
% model='itg-grace2010s'; %nmax 180
% model='GOCO03S';  %nmax 250
model='egm2008';  %nmax=2190

nmax=500;  %default grid step: 10 arcmin, see help to 'compute_geopot_grids'
% nmax=90;  %default grid step: 1 deg, see help to 'compute_geopot_grids'

% Computation of grid for the selected geopotential functional
[lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh');

%% Geoid height in 3D as PNG image
model='egm2008'; nmax=500;    %selection of geopotential model
[lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh');
[hc,hlab,name_png]=rotating_3d_globe(lond,latd,gh,'coastlines',1,...
   'exaggeration_factor',1.3e4,'radius',6378e3,'units','m',...
   'graph_label',sprintf('Geoid height (%s, nmax=%d)',upper(model),nmax),...
   'clbr_limits',[-90 90],'clbr_tick',-100:20:100,...
   'cptcmap_pm','BlueWhiteOrangeRed',...
   'preview_figure_visible',1,...
   'tn_width',200,...
   'window_height',650);

%% Geoid height in 3D as animated GIF: full revolution, smooth motion
model='egm2008'; nmax=180;    %selection of geopotential model
[lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh');
[hc,hlab,name_png]=rotating_3d_globe(lond,latd,gh,...
   'radius',6378e3,'units','m',...
   'exaggeration_factor',1.3e4,...
   'coastlines',1,...
   'cptcmap_pm','BlueWhiteOrangeRed',...
   'clbr_limits',[-90 90],'clbr_tick',-100:20:100,...
   'preview_figure_visible',0,...
   'clbr_anim',0,...
   'anim_gif',1,...
   'anim_angle',360,'time_for_360deg',42,'fps',7,...
   'window_height_intermediate',330,...
   'window_height',110);


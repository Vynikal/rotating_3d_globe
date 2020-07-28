% m-file: driver_rotating_3d_globe_9_Manual_03_anim_Earth_topography
% Examples in the Manual: 3. Making of 3D animations step by step
% Examples in the Manual: 4.3. Geoid height in 3D as animated GIF/compressed video: 1 fps
% folder: ../fig_large
%
% When using this driver or fuctions it calls, please cite the reference:
%  Bezdek A, Sebera J, 2013. MATLAB script for 3D visualizing geodata on a rotating globe.
%  Computers & geosciences 56, 127�130. http://dx.doi.org/10.1016/j.cageo.2013.03.007

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015
% 6/2016 update to new Matlab versions up to R2015a

%% Earth topography in 3D as animated GIF: 36-degree segment, 1 fps
model='ETOPO2_010arcmin'; load (model);
rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'graph_label',sprintf('Earth topography: (%s)',model),...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'anim_gif',1,...
   'anim_angle',36,'time_for_360deg',42,'fps',1,...
   'window_height',350);
%% Earth topography in 3D as animated GIF: 36-degree segment, 10 fps
model='ETOPO2_010arcmin'; load (model);
rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'graph_label',sprintf('Earth topography: (%s)',model),...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'anim_gif',1,...
   'anim_angle',36,'time_for_360deg',42,'fps',10,...
   'window_height',350);
%% Earth topography in 3D as animated GIF: 36-degree segment, 20 fps, flickering
model='ETOPO2_010arcmin'; load (model);
rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'graph_label',sprintf('Earth topography: (%s)',model),...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'anim_gif',1,...
   'anim_angle',36,'time_for_360deg',42,'fps',20,...
   'window_height',350);
%% Earth topography in 3D as animated GIF: 36-degree segment, 20 fps, smoothed
model='ETOPO2_010arcmin'; load (model);
rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'graph_label',sprintf('Earth topography: (%s)',model),...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'anim_gif',1,...
   'anim_angle',36,'time_for_360deg',42,'fps',20,...
   'window_height_intermediate',700,...
   'font_size',18,...
   'window_height',350);

%% Earth topography in 3D as animated GIF: full revolution, smooth motion
model='ETOPO2_030arcmin'; load (model);
[hc,hlab,filename_anim]=rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'preview_figure_visible',0,...
   'clbr_anim',0,...
   'anim_gif',1,...
   'anim_angle',360,'time_for_360deg',42,'fps',7,...
   'window_height_intermediate',330,...
   'window_height',110);

%% Earth topography in 3D animation: PNG files for slideshow
model='ETOPO2_010arcmin'; load (model);
rotating_3d_globe(lond_etopo2,latd_etopo2,elev_etopo2_km,...
   'exaggeration_factor',13,'radius',6378,'units','km',...
   'cptcmap_pm','GMT_globe',...
   'clbr_limits',[-10 10],'clbr_tick',-10:2:10,...
   'graph_label',sprintf('Earth topography: (%s)',model),...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'anim_gif',2,...
   'anim_angle',360,'time_for_360deg',36,'fps',1,...
   'window_height',500);


%% Geoid height in 3D as animated GIF/compressed video: 1 fps
% jako ze to je vhodne napr. do prezentace
model='egm2008'; nmax=500;    %selection of geopotential model
[lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh');
[hc,hlab,name_png]=rotating_3d_globe(lond,latd,gh,...
   'coastlines',1,...
   'exaggeration_factor',1.3e4,'radius',6378e3,'units','m',...
   'cptcmap_pm','BlueWhiteOrangeRed',...
   'graph_label',sprintf('Geoid height (%s, nmax=%d)',upper(model),nmax),...
   'clbr_limits',[-90 90],'clbr_tick',-100:20:100,...
   'preview_figure_visible',0,...
   'tn_width',200,...
   'clbr_anim',1,'clbr_reduce',0.4,...
   'video_format','wmv',...
   'anim_gif',1,...
   'anim_angle',360,'time_for_360deg',42,'fps',1,...
   'window_height',500);


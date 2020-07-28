% m-file: test_file_for_3d_printer
% A test driver script to create an STL file appropriate for 3D printers.
% You need the package 'rotating_3d_globe.zip' to be installed.
% It takes me a few seconds on my laptop to produce a smooth Earth (nmax=30;).
% To produce a realistic model, uncomment nmax=500 below, then producing the STL file took 6 minutes on my laptop.
% Best regards and good luck, Ales
% 1/2017, bezdek@asu.cas.cz

clear

%% Selection of geopotential model and computation/loading of the grid values
model='egm2008';  %nmax=2190
nmax=30;
% nmax=500;


% Computation of grid for the selected geopotential functional
[lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh');

%% Geoid height in 3D as PNG image
   [hc,hlab,name_png]=rotating_3d_globe(lond,latd,gh,'coastlines',1,...
      'exaggeration_factor',1.3e4,'radius',6378e3,'units','m',...
      'graph_label',sprintf('Geoid height (%s, nmax=%d)',upper(model),nmax),...
      'clbr_limits',[-90 90],'clbr_tick',-100:20:100,...
      'cptcmap_pm','BlueWhiteOrangeRed',...
      'preview_figure_visible',1,...
      'stl_file_export',1,...
      'window_height',450);

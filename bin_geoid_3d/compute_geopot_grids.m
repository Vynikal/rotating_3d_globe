function [lond,latd,geopot_grid]=compute_geopot_grids(model,nmax,varargin)
% [LOND,LATD,GEOPOT_GRID]=COMPUTE_GEOPOT_GRIDS(MODEL,NMAX,VARARGIN)
%     Computation/loading of grid for given geopotential functional using specified gravity field model
%
% 'model'..name of a mat file containing harmonic coefficients in matrices cnm and snm, which 
%     must contain coefficients of degree 'nmax' or higher
%     'model' is only used as a model label, if the specific spherical harmonic coefficients are input as optional arguments 'cnm' 'snm'
% 'nmax'..specify the maximum degree of the model to be used in the harmonic synthesis.
%     The higher 'nmax', the longer computation, but also the finer grid.
%
% Optional argument 'functional' can be one of the following: 
%    'gh'..geoid height (default); 
%    'gd'..gravity disturbance (negative radial first derivative); 
%    'ga'..gravity anomaly; 
%    'Trr'..second radial derivative
%    All the quantities are computed in the spherical approximation (Hofmann-Wellenhof Moritz 2005 Physical Geodesy, pp. 96-98).
%
% 'grid_stepd' step of the grid (deg)
%     The default value of 'grid_stepd' is rounded in arc-minutes to give 4 points per max. sine component, namely:
%          nmax=90     grid_stepd*60 = 60 arcmin = 1 deg
%          nmax=180    grid_stepd*60 = 30 arcmin = 0.5 deg
%          nmax=500    grid_stepd*60 = 10 arcmin = 0.1666... deg
%          nmax=1000   grid_stepd*60 =  5 arcmin = 0.8333... deg
%          nmax=1100   grid_stepd*60 =  4 arcmin = 0.8333... deg
%          nmax=2190   grid_stepd*60 =  2 arcmin = 0.0333... deg
% Specify 'grid_stepd' if you wish a particular value of grid step, e.g. for 2-deg grid of EGM2008 up to nmax=10:
%     [lond,latd,gh]=compute_geopot_grids('egm2008',10,'grid_stepd',2,'functional','gh');
% 
% 'subtract_normal_field'=0 .. no normal field subtracted
% 'subtract_normal_field'=1 .. the GRS 1980 normal field is subtracted (default) 
% 'subtract_normal_field'=2 .. the central term C00 is put equal to zero
% 'subtract_normal_field'=3 .. coefficients C00, C20 are put equal to zero
%
% 'cnm'..matrix of spherical harmonic coefficients
% 'snm'..matrix of spherical harmonic coefficients
% 'Re'=reference radius of the body (m)
%     6378136.3 m for Earth (default)/other value 
%
% 'gm'=geocentric gravitational constant (m3/s2)
%     3.986004415e14 m3/s2 for Earth (default)/other value 
%
% The function tests, whether the grid has already been computed and stored.
% In the folder 'temp', computed geopotential grids are stored for speeding up the execution of the function.
% You may delete its contents if you wish.
% We save only grids which took more than 'time_limit_to_save_grids' seconds to compute (default 1 sec).

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015
%% Default parameters
functional='gh';
grav_const
subtract_normal_field=1;
grid_stepd=floor(360/nmax/4*60)/60;
time_limit_to_save_grids=1;
%% Parsing optional arguments
i=1;
while i<=length(varargin)
   switch lower(varargin{i})
      case 'functional'
         functional=varargin{i+1};
      case 'grid_stepd'
         grid_stepd=varargin{i+1};
      case 'subtract_normal_field'
         subtract_normal_field=varargin{i+1};
      case 're'
         Re=varargin{i+1};
      case 'gm'
         gm=varargin{i+1};
      case 'cnm'
         cnm=varargin{i+1};
      case 'snm'
         snm=varargin{i+1};
   end;
   i=i+2;
end;
%% Check the path
% a simple way to check, whether also the subdirectories were added
if ~exist('asu-ch-0309.gfc','file')
   mfile= mfilename('fullpath');
   [pathstr, name, ext] = fileparts(mfile);
   addpath(genpath(pathstr));
end
%% Setting up 
grid_stepm=grid_stepd*60;  %arcmin
nmax1=nmax+1;
if ~exist('cnm','var')
   eval(['load ' model ]);
end

nmax_cnm=length(cnm)-1;
if nmax_cnm<nmax
   nmax=nmax_cnm;
   nmax1=nmax+1;
end

tic
fprintf('\n-------------------------------\n');
fprintf('Model: %s, nmax=%d, half-wave: %.3g km\n',model,nmax,20e3/nmax);
% fprintf('  half-wave: %.3g km\n',20e3/nmax);
fprintf('Grid step: %.3g deg = %.4g arcmin = %.3g km\n',grid_stepd,grid_stepm,grid_stepd/360*40e3);

if exist('functional','var')
   if strcmpi(functional,'gh')
      n_functional=1;
      label='geoid heights (m)';
   elseif strcmpi(functional,'ga')
      n_functional=2;
      label='gravity anomaly (m/s2)';
   elseif strcmpi(functional,'gd')
      n_functional=3;
      label='gravity disturbance (m/s2)';
   elseif strcmpi(functional,'Trr')
      n_functional=4;
      label='second radial derivative (1/s2)';
   end
end

fprintf('Compute_geopot_grids: %s\n',label);

if subtract_normal_field==1
   cnm(1:9,1:9)=cnm(1:9,1:9)-cnm_normal(8);
   fprintf('Normal field was subtracted.\n');
elseif subtract_normal_field==2
   cnm(1,1)=0;
   fprintf('The central term C00 was zeroed.\n');
elseif subtract_normal_field==3
   cnm(1,1)=0;
   cnm(3,1)=0;
   fprintf('Coefficients C00, C20 were put equal to zero.\n');
else
   fprintf('Normal field was not subtracted.\n');
end

Nlon=360/grid_stepd; Nlat=180/grid_stepd;

latd=90-(0:Nlat)*180/Nlat;
Nlat=Nlat+1;
colatd=90-latd;
Ngrid=Nlon*Nlat;
fprintf('Number of points: %d = %.1f mil = %.1f mld\n',Ngrid,Ngrid/1e6,Ngrid/1e9);

%% Test, whether the grid is already computed
% In the folder 'temp' the computed geopotential grids are stored for speeding up the execution of the function.
% You may delete its contents if you wish.
% We save only grids which took more than 'time_limit_to_save_grids' seconds to compute.

filename=mfilename('fullpath');
[pathstr, name, ext] = fileparts(filename);
temp_folder=[pathstr '/temp/'];
if ~exist(temp_folder,'file')
   mkdir(temp_folder);
   fid=fopen([temp_folder 'readme.txt'],'w');
   
   fprintf(fid,'%% This folder was created by function: compute_geopot_grids.m\n');
   fprintf(fid,'%% In the folder ''temp'' the computed geopotential grids are stored for speeding up the execution of the function.\n');
   fprintf(fid,'%% You may delete its contents if you wish.\n');
   fprintf(fid,'%% We save only grids which took more than ''time_limit_to_save_grids'' seconds to compute.\n');
   fclose(fid);
end
addpath(temp_folder);

filename=sprintf('grid_%s_%d_%s_%03.3gm_subtr%d.mat',model,nmax,functional,grid_stepm,subtract_normal_field);
if exist(filename,'file')
   load(filename);
   fprintf('Loaded file: %s\n',filename);
   return;
end

%% Computation
h_waitbar= waitbar(0,sprintf('Computation: %s ...',label));

Am=zeros(Nlat,nmax1);
Bm=zeros(Nlat,nmax1);

for m=0:nmax
   m1=m+1;
   Pnm=alfs_Pnm_column(nmax,m,colatd);
   if n_functional==1      %label='geoid heights (m)';
      Am(:,m1)=Pnm*cnm(m1:nmax1,m1);
      Bm(:,m1)=Pnm*snm(m1:nmax1,m1);
   elseif n_functional==2  %label='gravity anomaly (m/s2)';
      Am(:,m1)=Pnm*((m-1:nmax-1)'.*cnm(m1:nmax1,m1));
      Bm(:,m1)=Pnm*((m-1:nmax-1)'.*snm(m1:nmax1,m1));
   elseif n_functional==3  %label='gravity disturbance (m/s2)';
      Am(:,m1)=Pnm*((m+1:nmax+1)'.*cnm(m1:nmax1,m1));
      Bm(:,m1)=Pnm*((m+1:nmax+1)'.*snm(m1:nmax1,m1));
   elseif n_functional==4  %label='second radial derivative (1/s2)';
      Am(:,m1)=Pnm*((m+1:nmax+1)'.*(m1+2:nmax1+2)'.*cnm(m1:nmax1,m1));
      Bm(:,m1)=Pnm*((m+1:nmax+1)'.*(m1+2:nmax1+2)'.*snm(m1:nmax1,m1));
   end
   waitbar((m+1)/nmax);
end

lond=(0:Nlon)*360/Nlon;
lond=lond-180;
lon=lond/rad;
cos_mlon=cos((0:nmax)'*lon);
sin_mlon=sin((0:nmax)'*lon);
if n_functional==1
   geopot_grid=Re*(Am*cos_mlon+Bm*sin_mlon);
elseif n_functional==2 || n_functional==3
   geopot_grid=gm/Re^2*(Am*cos_mlon+Bm*sin_mlon);
elseif n_functional==4
   geopot_grid=gm/Re^3*(Am*cos_mlon+Bm*sin_mlon);
end
toc1=toc;
close(h_waitbar);
fprintf('Time of computation: %.0f sec = %.1f min\n',toc1,toc1/60);

if toc1>time_limit_to_save_grids
   save([temp_folder filename],'latd','lond','geopot_grid');
   fprintf('Saved file: %s\n',filename);
end


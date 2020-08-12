% m-file: srov_cas_01
% 01a: pridame k nazvu modelu epochu kvuli pocitani gridu (model_epocha)
% skript k otestovani nacitani modelu s casovou zavislosti
% funguje pro Eigen6s, cilem je fungovani pro Eigen6S2
%  pro Eigen6S2 nutno upravit nacitaci fci: cti_icgem_model

clear all

% ref_grav_model_name='eigen-6s';
ref_grav_model_name='EIGEN-6S4_v2';

nmax=20;
%% shc pro 10/2005
% epocha ve zlomcich roku, napr. leden 2005: epoch_yr=2005+31/2
y=2005; m=10; d=15;
epoch_yr=jd2yr(cal2jd(y,m,d));
dtn=datenum(y,m,d);
model_epocha=sprintf('%s_%s',upper(ref_grav_model_name),datestr(dtn,'mmm yyyy'));
ymd=y*10000+m*100+d;
[cnm, snm, header]=cti_icgem_model(lower(ref_grav_model_name),nmax,ymd,'period',0.5);
model=header.modelname;

%% mapka gh
if 0
   % kod podle http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/index.php#1b
   
   % Computation of grid for the selected geopotential functional
   [lond,latd,gh]=compute_geopot_grids(model,nmax,'functional','gh',...
      'cnm',cnm,'snm',snm,'re',header.radius,'gm',header.earth_gravity_constant);
   % 'cnm'..matrix of spherical harmonic coefficients
   % 'snm'..matrix of spherical harmonic coefficients
   % 'Re'=reference radius of the body (m)
   % 'gm'=geocentric gravitational constant (m3/s2)
   
   % Geoid height as 2D map
   [hc,htit,name_png]=elevation_2d_map(lond,latd,gh,'coastlines',1,...
      'units','m',...
      'graph_label',sprintf('Geoid height (%s, nmax=%d,%s)',upper(model),nmax,datestr(dtn,'mmm yyyy')),...
      'clbr_limits',[-90 90],'clbr_tick',-100:20:100,...
      'cptcmap_pm','BlueWhiteOrangeRed',...
      'window_height',650);
    return
end
%% mapka jen pro TVG
% y=2005; m=1; d=1; ymd=y*10000+m*100+d;
[cnm_noSeas, snm_noSeas]=cti_icgem_model(lower(ref_grav_model_name),nmax,ymd,'seasonal',0);
% odecteme slozku konst+trend, zbude nam sinusova tvg, rijnova tvg, ta je treba na obr. 8 v HP
cnm_tvg=cnm-cnm_noSeas;
snm_tvg=snm-snm_noSeas;

% Computation of grid for the selected geopotential functional
[lond,latd,gh]=compute_geopot_grids(model_epocha,nmax,'functional','gh',...
   'subtract_normal_field',0,...
   'cnm',cnm_tvg,'snm',snm_tvg,'re',header.radius,'gm',header.earth_gravity_constant);
% pribyl zde pozadavek, aby se neodecetlo normalni pole

% v milimetrech
gh_mm=gh*1000;

% Geoid height as 2D map
[hc,htit,name_png]=elevation_2d_map(lond,latd,gh_mm,'coastlines',1,...
   'units','mm',...
   'graph_label',sprintf('Sezonni maximum TVG (GH, %s, nmax=%d)',model_epocha,nmax),...
   'clbr_limits',[-9 9],'clbr_tick',-100:2:100,...
   'cptcmap_pm','BlueWhiteOrangeRed',...
   'window_height',650);

soub=sprintf('srov_cas_01a_%s_%d_%d_%d.png',model,nmax,y,m);
movefile(name_png,soub);
crop(soub)


function [hc,hlab,name_png]=rotating_3d_globe(lond,latd,elev,varargin)
% ROTATING_3D_GLOBE  Produces 3D image/animation of elevation data given on a longitude-latitude grid.
%
% [HC,HLAB,NAME_PNG]=ROTATING_3D_GLOBE(LOND,LATD,ELEV,VARARGIN)
% INPUT:
% 'lond' N×1 vector of longitudes, 'latd' M×1 vector of latitudes, 'elev' N×M matrix of elevation data
% OUTPUT:
% Every time, a preview 3D image with a colorbar is produced.
% To create file(s) with animation, the user has to enter the necessary input for the animation.
% 'hc'..handle to colorbar; 'hlab'..text label; 'name_png'..name of the preview png file.
%
% For explanation and examples of use, see
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/rotating_3d_globe/manual.html
%
%
% OPTIONAL ARGUMENTS
% 'radius'..defines the radius of 3D globe to which the 'elev' data are added for representing the relief.
% 'exaggeration_factor'..how much the 'elev' data will be exaggerated relative to 'radius' according
%    to the formula: r=radius/exaggeration_factor+elev
% 'units'..will be given as label for colorbar, e.g. 'km'
% 'graph_label'..descriptive string displayed in figure and also in the name of the output image/animation
%
% 'cptcmap_pm'..name of color palette file (*.cpt)
% 'clrmap'..name of color map file, either built-in Matlab, or provided by user
% 'clbr_limits'..colorbar limits (important for figures well structured in colours)
% 'clbr_reduce'..reduction of size of displayed colorbar (default 0.6)
% 'clbr_tick'..vector of colorbar ticks
% 'clbr_label_position'=0/1..below the colorbar (0, default), to the right (1, only for Matlab starting with R2014b)
%
% 'coastlines'=0/1..display coastlines (1) or not (0, default).
% 'coastlines_lw'..line width of displayed coastlines
% 'countries'=0/1..display country boundaries
% 'countries_lw'..line width of displayed country boundaries
%
% 'grid'=0/1..display grid lines on the globe
% 'grid_lw'..line width of displayed grid lines
% 'grid_lond'..grid line step in longitude (deg)(default 60 deg)
% 'grid_latd'..grid line step in latitude (deg) (default 30 deg)
%
% 'font_size'..user defined font size
% 'background_color'..background color, e.g. 'k'
% 'elev_colour'..if only color shades given by 'elev' data are to projected on a perfect sphere
%
% 'window_height'..user defined size of output images/animation (in pixels; default 500)
% 'azimuth'..initial azimuth of the viewer (deg; default 110 deg)
% 'elevation'..initial elevation of the viewer (deg; default 10 deg)
% 'view_angle'..zooming; try values: 3 corresponds to regional view, 1 or 2 gets you more closer
%
% 'anim_angle'..limit the rotation angle, which is animated (deg; full angle 360 deg, for debugging use e.g. 36 deg)
% 'fps'..frames per second (default 1)
% 'time_for_360deg'..time of animation to make a full rotation of 360 deg (in sec; default 42)
% 'anim_gif'=0/1/2: 1..animation saved as animated gif; 2..individual png files for slideshow saved in a folder
% 'video_format'='wmv','mp4','xvid' animation saved as a compressed video file;
% 'video_format'='uncomp' animation saved as uncompressed Matlab avi video file.
%
% 'clbr_prev'=0/1..show the colorbar in the preview image (default 1)
% 'clbr_anim'=0/1..show the colorbar in the animation (default 0)
% 'preview_figure_visible'=0/1..show the figure window of the preview image (default 1)
% 'anim_preview_figure_visible'=0/1..show the figure window of animation (default 0)
% 'window_height_intermediate'..for smooth rotation with high fps, use an intermediate larger window
%
% 'tn_width' generates a thumbnail of a given width (e.g. thumb_width'..200); no action for thumb_width'..0;
% 'rend'..explicitely defined renderer. Sometimes under linux creation of images based on large data sets
%    crashed with the segmentation fault error, then setting the renderer to 'painters' resolved the problem.
%
% 'stl_file_export'=0/1..generates an STL file for 3D printers, see the online manual for an example script (section 6.2)
%
% For examples of using all optional arguments, please see:
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/rotating_3d_globe/manual.html

% Ales Bezdek, bezdek@asu.cas.cz
% 19/9/2013 longitude/latitude grid lines added
% 25/2/2016 solved problems: hardcopy warning removed, animated gifs problem solved (cq not needed), avifile replaced by VideoWriter
% 25/2/2016 clbr_label_position: starting with R2014b, graphic objects changed
%% Default parameter values
exaggeration_factor=1;
radius=1;
units='m';
graph_label='';
coastlines=0;
coastlines_lw='';
countries=0;
countries_lw='';
gridd=0;
gridd_lw='';
gridd_lond=60;
gridd_latd=30;

cptcmap_pm=''; %default colormap is the Matlab one
clrmap='';  %default colormap is the Matlab one
background_color='w';  % background color
clbr_limits='auto'; %colorbar limits (important for figures well structured in colours)
clbr_reduce=0.6;
font_size=0;
clbr_tick='';
clbr_label_position=0;

window_height=500;
azimuth0=110; elevation0=10; %specify the (initial) viewpoint in terms of azimuth and elevation
azimuth='';elevation='';
view_angle='';

video_format=''; % video_format=''..no video file produced; other choices: 'wmv','mp4','xvid'
anim_gif=0;

clbr_prev='';
clbr_anim=0;  %colorbar in animation
anim_preview_figure_visible=0;
preview_figure_visible=1;

anim_angle=36;
fps=1;
time_for_360deg=42;
elev_colour='';
window_height_intermediate=0;

tn_width=0;  %generates a thumbnail of a given width (e.g. thumb_width=200); no action for thumb_width=0;

% Sometimes under linux creation of images based on large data sets crashed with the Segmentation fault error,
% then setting the renderer to 'painters' resolved the problem.
rend='';  %renderer
stl_file_export=0;
%% Parsing optional arguments
i=1;
while i<=length(varargin)
   switch lower(varargin{i})
      case 'exaggeration_factor'
         exaggeration_factor=varargin{i+1};
      case 'radius'
         radius=varargin{i+1};
      case 'units'
         units=varargin{i+1};
      case 'graph_label'
         graph_label=varargin{i+1};
      case 'coastlines'
         coastlines=varargin{i+1};
      case 'coastlines_lw'
         coastlines_lw=varargin{i+1};
      case 'countries'
         countries=varargin{i+1};
      case 'countries_lw'
         countries_lw=varargin{i+1};
      case 'grid'
         gridd=varargin{i+1};
      case 'grid_lw'
         gridd_lw=varargin{i+1};
      case 'grid_lond'
         gridd_lond=varargin{i+1};
      case 'grid_latd'
         gridd_latd=varargin{i+1};

      case 'cptcmap_pm'
         cptcmap_pm=varargin{i+1};
      case 'clrmap'
         clrmap=varargin{i+1};
      case 'background_color'
         background_color=varargin{i+1};

      case 'clbr_limits'
         clbr_limits=varargin{i+1};
      case 'clbr_tick'
         clbr_tick=varargin{i+1};
      case 'clbr_label_position'
         clbr_label_position=varargin{i+1};
      case 'clbr_reduce'
         clbr_reduce=varargin{i+1};
      case 'font_size'
         font_size=varargin{i+1};
      case 'elev_colour'
         elev_colour=varargin{i+1};

      case 'window_height'
         window_height=varargin{i+1};
      case 'azimuth'
         azimuth=varargin{i+1};
      case 'elevation'
         elevation=varargin{i+1};
      case 'view_angle'
         view_angle=varargin{i+1};

      case 'video_format'
         video_format=varargin{i+1};
      case 'anim_gif'
         anim_gif=varargin{i+1};
      case 'clbr_anim'
         clbr_anim=varargin{i+1};
      case 'anim_preview_figure_visible'
         anim_preview_figure_visible=varargin{i+1};
      case 'preview_figure_visible'
         preview_figure_visible=varargin{i+1};
      case 'clbr_prev'
         clbr_prev=varargin{i+1};
      case 'window_height_intermediate'
         window_height_intermediate=varargin{i+1};

      case 'anim_angle'
         anim_angle=varargin{i+1};
      case 'time_for_360deg'
         time_for_360deg=varargin{i+1};
      case 'fps'
         fps=varargin{i+1};

      case 'tn_width'
         tn_width=varargin{i+1};

      case 'rend'
         rend=varargin{i+1};

      case 'stl_file_export'
         stl_file_export=varargin{i+1};

   end;
   i=i+2;
end;
%% Various
window_height0=window_height;
if window_height_intermediate
   window_height=window_height_intermediate;
end
if isempty(azimuth)
   azimuth=azimuth0;
end
if isempty(elevation)
   elevation=elevation0;
end
hc=[]; hlab=[]; name_png=[];
%% Check the path
% a simple way to check, whether also the subdirectories were added
if ~exist('asu-ch-0309.gfc','file')
   mfile= mfilename('fullpath');
   [pathstr, name, ext] = fileparts(mfile);
   addpath(genpath(pathstr));
end
%% Computation of spherical data for 3D image
% spherical coordinates
rad=57.295779513082323;
theta=(90-latd)/rad;
phi=lond/rad;
% lond, latd may already be given as meshgrids
if ~(sum(size(phi)==size(theta))==2 && sum(size(phi)==size(elev))==2)
   [phi,theta]=meshgrid(phi,theta);
end

r=radius/exaggeration_factor+elev;
x=r.*sin(theta).*cos(phi);
y=r.*sin(theta).*sin(phi);
z=r.*cos(theta);

%% Computation of spherical coordinates for coastlines
if coastlines
   load NaturalEarth_ne_50m_coastline
   ii=isnan(latd_NaturalEarth_ne_50m_coastline); %indices where NaN's are located
   latd1=latd_NaturalEarth_ne_50m_coastline;
   lond1=lond_NaturalEarth_ne_50m_coastline;
   latd1(ii)=0;
   lond1(ii)=0;

   gh_coast=interp2(lond,latd,elev,lond1,latd1);

   gh_coast(ii)=nan;
   latd1(ii)=nan;
   lond1(ii)=nan;

   theta=(90-latd1)/rad;
   phi=lond1/rad;

   r=radius/exaggeration_factor+gh_coast;
   xx=r.*sin(theta).*cos(phi);
   yy=r.*sin(theta).*sin(phi);
   zz=r.*cos(theta);
end

%% Computation of spherical coordinates for countries
if countries
   load NaturalEarth_ne_50m_countries
   ii=isnan(latd_NaturalEarth_ne_50m_countries); %indices where NaN's are located
   latd1=latd_NaturalEarth_ne_50m_countries;
   lond1=lond_NaturalEarth_ne_50m_countries;
   latd1(ii)=0;
   lond1(ii)=0;

   gh_coast=interp2(lond,latd,elev,lond1,latd1);

   gh_coast(ii)=nan;
   latd1(ii)=nan;
   lond1(ii)=nan;

   theta=(90-latd1)/rad;
   phi=lond1/rad;

   r=radius/exaggeration_factor+gh_coast;
   xx1=r.*sin(theta).*cos(phi);
   yy1=r.*sin(theta).*sin(phi);
   zz1=r.*cos(theta);
end
%% Computation of spherical coordinates for grid lines
if gridd
   % direction of grid lines: south-north
   dlat=.1;
   [lond2,latd2]=meshgrid(-180:gridd_lond:180-gridd_lond,-90:dlat:90+dlat);
   lond1=lond2(:);
   latd1=latd2(:);
   ii=latd1==90+dlat;
   latd1(ii)=0;
   lond1(ii)=0;

   gh_coast=interp2(lond,latd,elev,lond1,latd1);

   gh_coast(ii)=nan;
   latd1(ii)=nan;
   lond1(ii)=nan;

   theta=(90-latd1)/rad;
   phi=lond1/rad;

   r=radius/exaggeration_factor+gh_coast;
%    r=r*1.05;  %sometimes the grid-line dots are not visible, I don't know why, so I included this arbitrary factor
   r=r*1.0000001;  %sometimes the grid-line dots are not visible, I don't know why, so I included this arbitrary factor
   xxg=r.*sin(theta).*cos(phi);
   yyg=r.*sin(theta).*sin(phi);
   zzg=r.*cos(theta);

   % direction of grid lines: west-east
   dlon=2*dlat;
   [latd2,lond2]=meshgrid(-90+gridd_latd:gridd_latd:90-gridd_latd,-180:dlon:180+dlon);
   lond1=lond2(:);
   latd1=latd2(:);
   ii=lond1==180+dlon;
   latd1(ii)=0;
   lond1(ii)=0;

   gh_coast=interp2(lond,latd,elev,lond1,latd1);

   gh_coast(ii)=nan;
   latd1(ii)=nan;
   lond1(ii)=nan;

   theta=(90-latd1)/rad;
   phi=lond1/rad;

   r=radius/exaggeration_factor+gh_coast;
   r=r*1.002;  %sometimes the grid-line dots are not visible, I don't know why, so I included this arbitrary factor
   xxgg=r.*sin(theta).*cos(phi);
   yygg=r.*sin(theta).*sin(phi);
   zzgg=r.*cos(theta);
end


%% Figure windows: one figure for the preview, possibly the second figure for the animation
% Preview figure is produced every time, the colorbar is shown
graph_label_filename1=strrep(graph_label,':','');
if azimuth~=azimuth0 || elevation~=elevation0
   graph_label_filename1=sprintf('%s_%g_%g',graph_label_filename1,roundn(azimuth-90,-1),roundn(elevation,-1));
end
graph_label_filename1=strrep(graph_label_filename1,' ','_');
graph_label_filename1=strrep(graph_label_filename1,'(','');
graph_label_filename1=strrep(graph_label_filename1,')','');
graph_label_filename1=strrep(graph_label_filename1,'=','');
graph_label_filename1=strrep(graph_label_filename1,';','-');
graph_label_filename1=strrep(graph_label_filename1,'+','-');
graph_label_filename1=strrep(graph_label_filename1,':','-');
graph_label_filename1=strrep(graph_label_filename1,',','');
graph_label_filename1=strrep(graph_label_filename1,'^','');
graph_label_filename1=strrep(graph_label_filename1,'\','');
graph_label_filename1=strrep(graph_label_filename1,'/','');
graph_label_filename1=strrep(graph_label_filename1,'"','');
if ~isempty(view_angle)
   graph_label_filename1=sprintf('%s_va%g',graph_label_filename1,view_angle);
end
anim_avi=~isempty(video_format);
if anim_avi || anim_gif
   n_figures=2;
   %    time_for_360deg=36*(anim_angle/360);
   dt=1/fps;
   nFrames=ceil(time_for_360deg/dt*anim_angle/360);
   step_anim=anim_angle/nFrames;
   if ~window_height_intermediate
      filename_anim=sprintf('rotating_3d_globe_%s_px%04d_angle%g_fps%g',graph_label_filename1,window_height0,anim_angle,fps);
   else
      filename_anim=sprintf('rotating_3d_globe_%s_px%04d_px%04d_angle%g_fps%g',graph_label_filename1,window_height0,window_height_intermediate,anim_angle,fps);
   end
   filename_anim=strrep(filename_anim,'.','');
else
   n_figures=1; %only the 3D preview figure is produced
end

%% loop over figure(s) generation
close all
for i=1:n_figures
   if i==1
      if preview_figure_visible
         h_prev=figure('visible', 'on'); %default for preview
      else
         h_prev=figure('visible', 'off');
      end
      clbr=1;
      if ~isempty(clbr_prev)
         clbr=clbr_prev;
      end
      h_actual=h_prev;
   elseif i==2
      clbr=clbr_anim;
%       close all;
      h_waitbar=waitbar(0,'Animation is under way...');
      if anim_preview_figure_visible
         h_anim=figure('visible', 'on');
      else
         h_anim=figure('visible', 'off'); %default for animation
      end
      h_actual=h_anim;
      if anim_avi
         if verLessThan('matlab','7.12.0') %R2011a
            avi_file=avifile([filename_anim '.avi'],'fps',fps,'compression', 'None');
         else
            writerObj=VideoWriter([filename_anim '.avi'],'Uncompressed AVI');
            writerObj.FrameRate = fps;
            open(writerObj);            
         end
      end
      if anim_gif % Preallocate movie structure.
         mov(1:nFrames) = struct('cdata', [],'colormap', []);
      end
      name_png=filename_anim;
   end
   if ~isempty(rend); set(h_actual,'Renderer',rend); end
   set(h_actual,'color',background_color);
   % Setting up the graphical window
   % reduction of margins
   jednotky=get(h_actual,'Units');
   set(h_actual,'Units','normal')
   if clbr
      % If you wanted to leave a little space for axis labels, you could substitute:
      set(gca,'Position',[.05 .05 .9 .9])
      %       set(gca,'Position',[.1 .1 .8 .8])
   else
      set(gca,'Position',[0 0 1 1])
   end
   set(h_actual,'Units',jednotky);

   w=window_height;
   if window_height<=400; fs=9;
   elseif window_height>500; fs=12;
   else fs=10;
   end
   if isunix  %by experience, on my Ubuntu linux, the fonts are smaller by approx 1-2 points
      fs=fs+1;
   end
   if font_size, fs=font_size; end
   %    fprintf('fs=%d\n',fs);
   %    get(h_actual,'Renderer')

   if w>=1080
      w=w/(4/3);
   end
   if clbr; hw=w; w=w*1.2;
   else hw=w;
   end
   hw=ceil(hw); w=ceil(w);
   if rem(w,2)==1
      w=w-1;
   end
   if rem(hw,2)==1
      hw=hw-1;
   end
%    set(h_actual, 'WindowStyle','normal','Position',[10,10,w,hw],'PaperPositionMode','auto', 'PaperUnits','points','PaperSize',[w,hw]);
   set(h_actual,'WindowStyle','normal','Position',[10,10,w,hw],'PaperPositionMode','auto','PaperUnits','points','PaperSize',[w,hw],'PaperPosition',[0,0,w,hw]);
   if hw<1080
      set(h_actual,'PaperPositionMode','auto'); % saved figure size matches the displayed figure size
   else
      set(h_actual,'PaperPositionMode','manual'); % saved figure size differs from the display
   end

   if isempty(elev_colour)
      elev_colour=elev;
   end

   %% 3D rendering of elevation data
   surf(x,y,z,elev_colour,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
   if stl_file_export
      name_stl=sprintf('rotating_3d_globe_%s.stl',graph_label_filename1);
      fprintf('\nSTL file has been created:\n   %s\n',name_stl);
      surf2stl(name_stl,x,y,z);
   end
   camlight('right')
   material dull  %no specular reflections

   if ~isempty(cptcmap_pm)
      cmap = cptcmap(cptcmap_pm);
      colormap(cmap);
      graph_label_filename1=sprintf('%s_%s',graph_label_filename1,cptcmap_pm);
      if i==1 && n_figures==2
         filename_anim=sprintf('%s_%s',filename_anim,cptcmap_pm);
      end
   end

   if ~isempty(clrmap)
      colormap(clrmap);
      graph_label_filename1=sprintf('%s_%s',graph_label_filename1,clrmap);
      if i==1 && n_figures==2
         filename_anim=sprintf('%s_%s',filename_anim,clrmap);
      end
   end

   if clbr
      hc=colorbar('fontSize',fs);
      xlabel(hc,units,'fontSize',fs)
   end
   caxis(clbr_limits);
   if clbr && ~isempty(clbr_tick)
      set(hc,'ytick',clbr_tick);
   end

   if clbr && clbr_reduce
      pos=get(hc,'position');
      y1=pos(2)+(1-clbr_reduce)*pos(4)/2;
      if hw>1080
         set(hc,'position',[pos(1) y1 pos(3) pos(4)*clbr_reduce]);
      else
         set(hc,'position',[pos(1) y1 pos(3)*clbr_reduce pos(4)*clbr_reduce]);
      end
   end
  
   % I prefer to have the label below the colorbar as it was before R2014b
   if clbr && clbr_label_position==0 && ~verLessThan('matlab','8.4.0') %R2014b
      ylim1=hc.Limits;
      hc.Label.Rotation=0;   
      hc.Label.Position=[0.5 ylim1(1)-.03*abs(ylim1(1)) 0];
%       hc.Label.HorizontalAlignment='left';
   end
   
   %this is here by experimenting to have the Moon, which has a different
   %view, of the same size in the picture as other bodies
   view(110,10);

   axis equal
   axis vis3d
   axis off
   view(azimuth,elevation); %not to change the situation with respect to other views

   if coastlines
      hold on
      if window_height>4000
         lw=3;
      elseif window_height>2000
         lw=2;
      else
         lw=1;
      end
      if ~verLessThan('matlab','8.4.0') %R2014b
         lw=lw/2;    %new graphics engine changed the look of the graphs
      end      
      if ~isempty(coastlines_lw)
         lw=coastlines_lw;
      end
      plot3(xx,yy,zz,'k','linewidth',lw)
   end

   if countries
      hold on
      if window_height>4000
         lw=3;
      elseif window_height>2000
         lw=2;
      else
         lw=1;
      end
      if ~verLessThan('matlab','8.4.0') %R2014b
         lw=lw/2;    %new graphics engine changed the look of the graphs
      end      
      if ~isempty(countries_lw)
         lw=countries_lw;
      end
      plot3(xx1,yy1,zz1,'k','linewidth',lw)
   end

   if gridd
      hold on
      if window_height>4000
         lw=3;
      elseif window_height>2000
         lw=2;
      else
         lw=1;
      end
      if ~verLessThan('matlab','8.4.0') %R2014b
         lw=lw/2;    %new graphics engine changed the look of the graphs
      end      
      if ~isempty(gridd_lw)
         lw=gridd_lw;
      end
      plot3(xxg,yyg,zzg,'k:','linewidth',lw)
      plot3(xxgg,yygg,zzgg,'k:','linewidth',lw)
   end

   if ~isempty(view_angle)
      camva(view_angle);
   end

   h_light = findobj(h_actual,'Type','light');
   camlight(h_light,'right')

   graph_label_label=strrep(graph_label,'_','\_');

   if 1
      hlab=annotation(h_actual,'textbox',[.03 0.03 0.41 0.04667],...
         'FitBoxToText','on','LineStyle','none',...
         'String',{graph_label_label},...
         'fontSize',fs);
   else
      % to put the title above the figure
%      hlab=annotation(h_actual,'textbox',[.2 0.88 0.41 0.04667],...
      hlab=annotation(h_actual,'textbox',[0 0.88 1 0.04667],...
         'HorizontalAlignment','center',...
         'FitBoxToText','on','LineStyle','none',...
         'String',{graph_label_label},...
         'fontSize',fs);
   end
   
   % to preserve the background color in the saved figure
   set(h_actual,'InvertHardcopy','off');

   if i==1  %print out the preview
      fprintf('\nFunction: ROTATING_3D_GLOBE\n');
      fprintf('Elevation data:  min=%.3g, max=%.3g\n',min(min(elev)),max(max(elev)));
      if ~window_height_intermediate
         name_png=sprintf('rotating_3d_globe_preview_%s_px%04d.png',graph_label_filename1,window_height);
         eval(sprintf('print -dpng -r0 %s',name_png));
      else
         img=zbuffer_cdata(h_actual);
         zoom_factor=window_height0/window_height_intermediate;
         img_new=resize_image_ab(img,zoom_factor);
         name_png=sprintf('rotating_3d_globe_preview_%s_px%04d_px%04d.png',graph_label_filename1,window_height0,window_height_intermediate);
         imwrite(img_new,name_png);
      end
      fprintf('Preview PNG image has been created:\n   %s\n',name_png);
      %% Thumbnail generation
      if tn_width
         img=zbuffer_cdata(h_actual);
         zoom_factor1=tn_width/size(img,2);
%          zoom_factor1=tn_width/w;
         img_new=resize_image_ab(img,zoom_factor1);
         if ~verLessThan('matlab','8.4.0') %R2014b
            img_new(end,:,:)=[];
         end
         imwrite(img_new,sprintf('tn%d_%s.png',tn_width,name_png(1:end-4)),'png');
         fprintf('Thumbnail created:\n   tn%d_%s.png\n',tn_width,name_png(1:end-4));
      end
   end
   if i==2 %animation
      camlight(h_light,'right')
      for j=0:nFrames-1
         alpha=j*step_anim;
         %          fprintf('Frame generated at angle: %7.3f deg\n',alpha);
         camlight(h_light,'right')
         drawnow;
         if ~window_height_intermediate
            img=zbuffer_cdata(h_anim);
         else
            img1=zbuffer_cdata(h_anim);
            img=resize_image_ab(img1,zoom_factor);
         end
         if anim_avi
            if verLessThan('matlab','7.12.0') %R2011a
               F = im2frame(img);
               avi_file = addframe(avi_file, F);
            else
               writeVideo(writerObj,img);
            end
         end

         if anim_gif
            F = im2frame(img);
            mov(j+1)=F;
         end
         camorbit(-step_anim,0)
%          drawnow;
         waitbar(alpha/anim_angle,h_waitbar,sprintf('Animation is under way: %.0f %%',alpha/anim_angle*100));
      end
      close(h_waitbar);
   end
end
%% Conversion of uncompressed avi into the chosen video format
% for ffmpeg it is important that width and height of input video are even numbers, therefore we crop it if it is necessary
% but some binary distributions of ffmpeg does not recognize this crop syntax
if anim_avi
   if verLessThan('matlab','7.12.0') %R2011a
      avi_file=close(avi_file); %this is not mistake, it must be written this way
   else
      close(writerObj);
   end
   avi_error=...
      ['\nWARNING: A problem occured using external program "ffmpeg", the compressed video file cannot be produced,\n'...
      'only the huge uncompressed AVI video file has been created. \n'...
      'Consider installation of this free conversion tool, http://ffmpeg.org/.\n'...
      'After the installation, add its directory to the matlab path.\n'];

   if strcmpi(video_format,'wmv')
      filename_anim_ext=[filename_anim '.' video_format];
      if exist(filename_anim_ext,'file'); delete(filename_anim_ext); end
%       [status,result]=system(sprintf('ffmpeg -i %s.avi -vf crop=%d:%d:0:0 -qscale 1 %s.wmv',filename_anim,w,hw,filename_anim));
      [status,result]=system(sprintf('ffmpeg -i %s.avi -qscale 1 %s.wmv',filename_anim,filename_anim));
      if ~status
         delete([filename_anim '.avi']);
         fprintf('Compressed %s video file is created:\n   %s\n',upper(video_format),filename_anim_ext);
      else
         fprintf(avi_error);
         fprintf('SYSTEM ERROR MESSAGE: %s\n',result);
      end
   elseif strcmpi(video_format,'mp4')
      filename_anim_ext=[filename_anim '.' video_format];
      if exist(filename_anim_ext,'file'); delete(filename_anim_ext); end
%       [status,result]=system(sprintf('ffmpeg -i %s.avi -vf crop=%d:%d:0:0 -vcodec mpeg4 -qscale 1 %s.mp4',filename_anim,w,hw,filename_anim));
      [status,result]=system(sprintf('ffmpeg -i %s.avi -vcodec mpeg4 -qscale 1 %s.mp4',filename_anim,filename_anim));
      if ~status
         delete([filename_anim '.avi']);
         fprintf('Compressed %s video file is created:\n   %s\n',upper(video_format),filename_anim_ext);
      else
         fprintf(avi_error);
         fprintf('SYSTEM ERROR MESSAGE: %s\n',result);
      end
   elseif strcmpi(video_format,'xvid')
      filename_anim_ext=[filename_anim '.avi'];
      uiwait(h_actual,1);
      movefile([filename_anim '.avi'],[filename_anim '1.avi']);
%       [status,result]=system(sprintf('ffmpeg -i %s1.avi -vf crop=%d:%d:0:0 -vcodec mpeg4 -vtag xvid -qscale 1 %s.avi',filename_anim,w,hw,filename_anim));
      [status,result]=system(sprintf('ffmpeg -i %s1.avi -vcodec mpeg4 -vtag xvid -qscale 1 %s.avi',filename_anim,filename_anim));
      if ~status
         delete([filename_anim '1.avi']);
         fprintf('Compressed %s video file is created:\n   %s\n',upper(video_format),filename_anim_ext);
      else
         fprintf(avi_error);
         fprintf('SYSTEM ERROR MESSAGE: %s\n',result);
      end
   else
      filename_anim_ext=[filename_anim '.avi'];
      fprintf('Uncompressed AVI video file is created:\n   %s\n',filename_anim_ext);
   end
end
%% Animated gif
if anim_gif==2  %only the folder with individual png files
   if exist(filename_anim,'file'), delete([ filename_anim './*']); end
   if ~exist(filename_anim,'file'), mkdir(filename_anim); end
   for i=1:length(mov)
      img = frame2im(mov(i));
      eval(sprintf('imwrite(img,''%s/png_%06d.png'',''png'');',filename_anim,i));
   end
elseif anim_gif  %animated gif
   comp=computer;
   
   % Starting with R2011, the function rgb2ind should be included in the main Matlab module
   % so it is not necessary to include it separately
   if ~verLessThan('matlab','7.12.0') %R2011a
      x = path; % string
      y=strread(path,'%s','delimiter', pathsep);
      for i=1:length(y)
         if ~isempty(strfind(y{i},'graphics_functions_for_R2008a'))
%             fprintf('Removed from Matlab path: %s\n',y{i});
            rmpath(y{i});
         end
      end
   end
   
   % manually I copied the necessary mex functions for 64bit Win7 and Linux
%    if strcmpi(comp,'PCWIN64') || strcmpi(comp,'GLNXA64') || exist('rgb2ind','file')
   if exist('rgb2ind','file')
      fprintf('Animated GIF is produced:\n   %s\n',[filename_anim '.gif']);
      for i=1:length(mov)
         img = frame2im(mov(i));
         %       [imind,cm] = rgb2ind(im,256,'nodither');
         [imind,cm] = rgb2ind(img,256,'dither');
         if i == 1;
            imwrite(imind,cm,[filename_anim '.gif'],'gif','DelayTime',dt,'Loopcount',inf);
         else
            imwrite(imind,cm,[filename_anim '.gif'],'gif','DelayTime',dt,'WriteMode','append');
         end
      end
   else
      if exist(filename_anim,'file')
         delete(filename_anim);
      end
      mkdir(filename_anim);
      anim_gif_error=[...
         'WARNING: Animated GIF cannot be produced because the function "rgb2ind.m" is not on your Matlab path. \n'...
         'There are several options of what to do in order to obtain animated gif.\n'...
         sprintf('(1) Folder named "%s"\n',filename_anim)...
         'was created with individual PNG files, which you can convert into an animated GIF using an external tool.\n'...
         '(2) For older versions of Matlab (e.g. R2008a), you can install Image Processing Toolbox, \n'...
         'which contains the function "rgb2ind.m".\n'...
         '(3) For older versions of Matlab (e.g. R2008a), you can write to the author of this function (bezdek@asu.cas.cz),\n'...
         'and he might send you the required few functions provided he has them for your operating system.\n'...
         '(4) In case of newer version (e.g. R2011a), there should be no problem, the function "rgb2ind.m" \n'...
         'is part of the basic Matlab module. Still you can contact me, bezdek@asu.cas.cz.\n'];
      fprintf(anim_gif_error);
      for i=1:length(mov)
         im = frame2im(mov(i));
         eval(sprintf('imwrite(im,''%s/png_%06d.png'',''png'');',filename_anim,i));
      end
   end
end
%% Undocumented hardcopy function
function cdata = zbuffer_cdata(hfig)
% http://www.mathworks.com/support/solutions/en/data/1-3NMHJ5/?solution=1-3NMHJ5
% Get CDATA from hardcopy using zbuffer
% 25/2/2016 modified by Ales Bezdek to comply with Matlab versions starting with R2014b

% Need to have PaperPositionMode be auto
orig_mode = get(hfig, 'PaperPositionMode');
set(hfig, 'PaperPositionMode', 'auto');

if verLessThan('matlab','8.4.0') %R2014b
   fig_Renderer     = get(hfig, 'Renderer');
   % Set image driver:
   if strcmpi(fig_Renderer, 'painters')
      imageDriver = '-dzbuffer';
   else
      imageDriver = ['-d', fig_Renderer];
   end
   % cdata = hardcopy(hfig, '-Dzbuffer', '-r0');
   cdata = hardcopy(hfig,imageDriver, '-r0');
else
   % http://www.mathworks.com/matlabcentral/answers/99925-why-does-the-screensaver-get-captured-when-i-use-getframe-on-a-matlab-figure-window-while-creating-a
   cdata = print(hfig,'-RGBImage','-r0');
end

% Restore figure to original state
set(hfig, 'PaperPositionMode', orig_mode); % end

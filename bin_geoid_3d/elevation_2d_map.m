function [hc,htit,name_png]=elevation_2d_map(lond,latd,elev,varargin)
%ELEVATION_2D_MAP Two-dimensional map of elevations.
%
% [HC,HTIT,NAME_PNG]=ELEVATION_2D_MAP(LOND,LATD,ELEV,VARARGIN) draws a color map of M×N matrix of elevations elev
%     as function of N×1 vector of longitudes lond (degrees) and M×1 vector of latitudes (degrees).
%     It returns the handle 'hc' to the colorbar, handle 'htit' to the title,
%     'name_png' string with the png image filename
%
% For examples of use, please see
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%  http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/rotating_3d_globe/manual.html
%
% OPTIONAL ARGUMENTS
% 'units'..will be given as label for colorbar, e.g. 'km'
% 'graph_label'..title text displayed in figure and also in the name of the output image
%
% 'map_center'..centres the view point either on 'Europe' (default), on 'Pacific', on a given meridian
%     (give the longitude as a positive number in degrees, e.g. 'meridian_180', 'meridian_270.75', etc.)
%
% 'projection'..'Robinson', 'Mollweide' (default, because it is equal-area and thus better for presenting scientific data)
% 'cptcmap_pm'..name of color palette file (*.cpt)
% 'clrmap'..name of color map file, either built-in Matlab, or provided by user
% 'clbr_limits'..colorbar limits (important for figures well structured in colours)
% 'clbr_tick'..vector of colorbar ticks
% 'clbr_reduce'..reduction of size of displayed colorbar (default 0.8)
%
% 'coastlines'=0/1..display coastlines (1) or not (0, default).
% 'countries'=0/1..display country boundaries
%
% 'font_size'..user defined font size
% 'show_ylabels'=0/1..show latitude labels on y axis (0, default)
%
% 'window_height'..user defined size of output images/animation (in pixels; default 500)
% 'preview_figure_visible'=0/1..show the figure window of the preview image (default 1)
%
% 'tn_width' generates a thumbnail of a given width (e.g. thumb_width'..200); no action for thumb_width'..0;
% 'rend'..explicitely defined renderer. Sometimes under linux creation of images based on large data sets
%    crashed with the segmentation fault error, then setting the renderer to 'painters' resolved the problem.
%
% Use is made of excellent M_Map: A mapping package for Matlab, http://www2.ocgy.ubc.ca/~rich/map.html.

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015
% 25/2/2016 modified thumbnails creation
%% Default parameter values
units='m';
graph_label='';
coastlines=0;
countries=0;
map_center='Europe';
% map_center='Pacific';

cptcmap_pm=''; %default colormap is the Matlab one
clrmap='';  %default colormap is the Matlab one
clbr_limits='auto'; %colorbar limits (important for figures well structured in colours)
font_size=0;
show_ylabels=0;
clbr_reduce=.8;
clbr_tick='';
pprojection='Mollweide';

window_height=500;

preview_figure_visible=1;

tn_width=0;  %generates a thumbnail of a given width (e.g. thumb_width=200); no action for thumb_width=0;

% Sometimes under linux creation of images based on large data sets crashed with the Segmentation fault error,
% then setting the renderer to 'painters' resolved the problem.
rend='';  %renderer
%% Parsing optional arguments

i=1;
while i<=length(varargin),
   switch lower(varargin{i})
      case 'units'
         units=varargin{i+1};
      case 'graph_label'
         graph_label=varargin{i+1};
      case 'coastlines'
         coastlines=varargin{i+1};
      case 'countries'
         countries=varargin{i+1};
         
      case 'map_center'
         map_center=varargin{i+1};
         
      case 'cptcmap_pm'
         cptcmap_pm=varargin{i+1};
      case 'clrmap'
         clrmap=varargin{i+1};
         
      case 'clbr_limits'
         clbr_limits=varargin{i+1};
      case 'clbr_tick'
         clbr_tick=varargin{i+1};
      case 'projection'
         pprojection=varargin{i+1};
      case 'font_size'
         font_size=varargin{i+1};
      case 'show_ylabels'
         show_ylabels=varargin{i+1};
         
      case 'window_height'
         window_height=varargin{i+1};
         
      case 'preview_figure_visible'
         preview_figure_visible=varargin{i+1};
         
      case 'tn_width'
         tn_width=varargin{i+1};
         
      case 'rend'
         rend=varargin{i+1};
         
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
%%
close all
if preview_figure_visible
   h_actual=figure('visible', 'on');
else
   h_actual=figure('visible', 'off'); %default, much faster creation of png for large data sets
end
clbr=1;
if ~isempty(rend); set(h_actual,'Renderer',rend); end

% Setting up the graphical window
% reduction of margins
if 0
   jednotky=get(h_actual,'Units');
   set(h_actual,'Units','normal')
   if clbr
      % If you wanted to leave a little space for axis labels, you could substitute:
      set(gca,'Position',[.05 .05 .9 .9])
      %           set(gca,'Position',[.1 .1 .8 .8])
   else
      set(gca,'Position',[0 0 1 1])
   end
   set(h_actual,'Units',jednotky);
end

w=window_height;
if window_height<=400; fs=9;
elseif window_height>=500; fs=12;
else fs=10;
end
fs=fs+1;
if ~verLessThan('matlab','8.4.0') %R2014b
   fs=fs+3;
end
if isunix  %by experience, on my Ubuntu linux, the fonts are smaller by approx 1-2 points
   fs=fs+2;
end
% fprintf('fs=%d\n',fs);
if font_size, fs=font_size; end
if clbr; hw=w; w=round(w*2);
else hw=w;
end
set(h_actual, 'WindowStyle','normal','Position',[10,10,w,hw],'PaperPositionMode','auto', 'PaperUnits','points','PaperSize',[w,hw]);

% 2D rendering of elevation data
if strcmpi(map_center,'Pacific') || strcmpi(map_center,'America')
   clond=-150;
   m_proj(pprojection,'clongitude',clond);
   m_pcolor(lond,latd,elev)
   hold on
   m_pcolor(lond-360,latd,elev)
elseif strcmpi(map_center,'Europe')
   clond=0;
   m_proj(pprojection,'clongitude',clond);
   m_pcolor(lond,latd,elev)
elseif ~isempty(map_center) && strcmpi(map_center(1:8),'meridian')
   clond=str2double(map_center(10:end));
   m_proj(pprojection,'clongitude',clond);
   m_pcolor(lond,latd,elev)
   hold on
   m_pcolor(lond+360,latd,elev)
else
   m_proj(pprojection);
   m_pcolor(lond,latd,elev)
end
shading interp;

graph_label_filename1=strrep(graph_label,':','-');
graph_label_filename1=strrep(graph_label_filename1,' ','_');
graph_label_filename1=strrep(graph_label_filename1,'(','');
graph_label_filename1=strrep(graph_label_filename1,')','');
graph_label_filename1=strrep(graph_label_filename1,'=','');
graph_label_filename1=strrep(graph_label_filename1,',','');
graph_label_filename1=strrep(graph_label_filename1,';','-');

if ~isempty(cptcmap_pm)
   cmap = cptcmap(cptcmap_pm);
   colormap(cmap);
   graph_label_filename1=sprintf('%s_%s',graph_label_filename1,cptcmap_pm);
end

if ~isempty(clrmap)
   colormap(clrmap);
   graph_label_filename1=sprintf('%s_%s',graph_label_filename1,clrmap);
end

if countries
   load NaturalEarth_ne_50m_countries
   %    ii=isnan(latd_NaturalEarth_ne_50m_countries); %indices where NaN's are located
   hold on
   m_plot(lond_NaturalEarth_ne_50m_countries,latd_NaturalEarth_ne_50m_countries,'k-');
   coastlines=0;
end

if coastlines
   %       hold on
   m_coast('color','black');
end

% m_grid('xticklabels',[],'yticklabels',[],'ytick',[-90 -60 -30 0 30 60 90]) %no y labels
m_grid('ytick',3,'yticklabels',[],'xticklabels',[])   %no labels
m_grid('xtick',3,'yticklabels',[],'xticklabels',[])   %no labels

if show_ylabels
   m_grid('xaxis','middle','ytick',[-60 -30 30 60])
end

if clbr
   hc=colorbar('fontSize',fs);
   xlabel(hc,units,'fontSize',fs)
   caxis(clbr_limits);
end

if clbr && ~isempty(clbr_tick)
   set(hc,'ytick',clbr_tick);
end

if clbr && ~verLessThan('matlab','8.4.0') %R2014b
   % I prefer to have the label below the colorbar as it was before R2014b
   ylim1=hc.Limits;
   hc.Label.Rotation=0;
   hc.Label.Position=[0.5 ylim1(1)-.03*abs(ylim1(1)) 0];
end

if clbr && clbr_reduce
   % Store the current positions of the axes and the colorbar.
   axpos=get(gca,'position');
   cpos = get(hc,'position');
   
   % Change the colorbar width to half of its original width by adjusting the third element in cpos. Set the colorbar Position property to the updated position vector. Then, reset the axes to its original position so that it does not overlap the colorbar.
   cpos(2)= cpos(2)+(1-clbr_reduce)*cpos(4)/2;
   cpos(3) = clbr_reduce*cpos(3);
   cpos(4) = clbr_reduce*cpos(4);
   set(hc,'position',cpos)
   set(gca,'position',axpos)
end

% m_grid('xticklabels',[]);   %equator with no labels
% m_grid('ytick',[-60 -30 30 60]) %with y labels
%    if isempty(strfind(str,'no_labels'))
% m_grid('ytick',[-60 -30 30 60]) %with y labels
%    else
%       m_grid('xtick',3,'yticklabels',[],'xticklabels',[])   %equator with no labels
%    end

graph_label_label=strrep(graph_label,'_','\_');
htit=title(graph_label_label,'fontSize',fs);

if strcmpi(map_center,'Pacific') || strcmpi(map_center,'America')
   name_png=sprintf('elevation_2d_map_%s_px%04d_pacific.png',graph_label_filename1,window_height);
elseif strcmpi(map_center,'Europe')
   name_png=sprintf('elevation_2d_map_%s_px%04d.png',graph_label_filename1,window_height);
elseif ~isempty(map_center) && strcmpi(map_center(1:8),'meridian')
   name_png=sprintf('elevation_2d_map_%s_px%04d_%03g.png',graph_label_filename1,window_height,clond);
else
   name_png=sprintf('elevation_2d_map_%s_px%04d.png',graph_label_filename1,window_height);
end

eval(sprintf('print -dpng -r0 %s',name_png));
fprintf('\nFunction: ELEVATION_2D_MAP\n');
fprintf('Elevation data:  min=%.3g, max=%.3g\n',min(min(elev)),max(max(elev)));
fprintf('2D map as a PNG image has been created:\n   %s.png\n',name_png);

%% Thumbnail generation
if tn_width
   img=zbuffer_cdata(h_actual);
   zoom_factor1=tn_width/w;
   img_new=resize_image_ab(img,zoom_factor1);
   if ~verLessThan('matlab','8.4.0') %R2014b
      img_new(end,:,:)=[];
   end
   imwrite(img_new,sprintf('tn%d_%s.png',tn_width,name_png(1:end-4)),'png');
   fprintf('Thumbnail created:\n   tn%d_%s.png\n',tn_width,name_png(1:end-4));
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

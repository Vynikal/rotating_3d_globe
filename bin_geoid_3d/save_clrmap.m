function save_clrmap(name)
%SAVE_CLRMAP   Saves the colormap of the current figure to a file.
%
% SAVE_CLRMAP(NAME) saves the colormap to a file 'name.m' for easy later use.
% 
%  Use the saved colormap in the same way as the usual MATLAB colormaps. 
%  One can easily modify existing colorbars, or create new colorbars, using the 'colormapeditor'.
%  For example, to save a modified or newly created colormap, use the command
%
%        save_clrmap('foo');
%
%  To use your saved colormap 'foo' in another MATLAB session, to set the colormap of the current figure, type
%
%        colormap(foo)
%
%
% See also colormapeditor

% Ales Bezdek, bezdek@asu.cas.cz, 11/2012

if ~exist('name','var')
   error('Error: Please provide a filename for the saved colormap: save_clrmap NAME');
end

% this does not work in R2016b
%map=get(gcf,'Colormap');

map=colormap(gca);

fid=fopen([name '.m'],'w');


fprintf(fid,'function map=%s\n',name);
fprintf(fid,'%%%s   User defined colormap saved by command ''save_clrmap''\n',upper(name));
fprintf(fid,'%%\n');
fprintf(fid,'%%   To reset the colormap of the current figure, use it in the same way as usual Matlab colormaps:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%             colormap(%s)\n',name);
fprintf(fid,'%%\n');
fprintf(fid,'%%   See also SAVE_CLRMAP\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%This help was automatically generated by function SAVE_CLRMAP.\n');
fprintf(fid,'\n');
% fprintf(fid,'%%Ales Bezdek, bezdek@asu.cas.cz, 11/2012\n');
% fprintf(fid,'%%\n');
% fprintf(fid,'%%\n');
fprintf(fid,'map=[\n');
for i=1:size(map,1)
   fprintf(fid,'%24.15e %24.15e %24.15e;\n',map(i,:));
end
fprintf(fid,'];\n');
fclose(fid);
function dtn=doy2dtn(yr,doy)
% dtn=doy2dtn(yr,doy)
% dtn..Matlab datenumber, yr..year, doy..day of year

% 3/2015 bezdek@asu.cas.cz
% according to: http://www.mathworks.com/matlabcentral/newsreader/view_thread/32006

dtn = datenum(yr-1, 12, 31, 0, 0, 0) + doy;


function dtn=ymd2dtn(ymd)
% dtn=ymd2dtn(ymd) converts YYYYMMDD format into Matlab serial date number
%
% Back: ymd=datestr(dtn,'yyyymmdd');
eval(sprintf('dtn=datenum(%s,%s,%s);',ymd(1:4),ymd(5:6),ymd(7:8)));


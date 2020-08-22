function dtn=ymd2dtn(ymd)
% dtn=ymd2dtn(ymd) converts YYYYMMDD.DDDD format into Matlab serial date number
%
% Back: ymd=datestr(dtn,'yyyymmdd');
str = num2str(ymd);
eval(sprintf('dtn=datenum(%s,%s,%s);',str(1:4),str(5:6),str(7:end)));


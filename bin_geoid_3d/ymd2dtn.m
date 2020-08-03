function dtn=ymd2dtn(ymd)
% dtn=ymd2dtn(ymd) converts YYYY_MMDD format into Matlab serial date number
%
% Back: ymd=datestr(dtn,'yyyy_mmdd');
eval(sprintf('dtn=datenum(%s,%s,%s);',ymd(1:4),ymd(6:7),ymd(8:9)));


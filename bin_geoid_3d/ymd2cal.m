function [y,m,d]=ymd2cal(ymd)
% [yr,mn,dy]=ymd2cal(ymd) converts YYYY.MMDD format into Y,M,D
% 
y=floor(ymd); 
m=floor((ymd-y)*1e2); 
d=(((ymd-y)*1e2)-m)*1e2;

function yr=jd2yr(jd)
% JD2YR  Converts Julian date to year and decimal of year.
% See also CAL2JD, DOY2JD, GPS2JD, JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR.
% Usage:   yr=jd2yr(jd)
% Input:   jd - Julian date
% Output:  yr - year and decimal of year

% Version: 24 Apr 99
% Simply vectorized by Ales Bezdek: 10/6/2008
if nargin ~= 1
  warning('Incorrect number of arguments');
  return;
end
if jd < 0
  warning('Julian date must be greater than or equal to zero');
  return;
end

[iyr,mn,yr] = jd2cal(jd);
jd0 = cal2jd(iyr,1,1);
jd1 = cal2jd(iyr+1,1,1);
%ab yr = iyr + (jd-jd0)/(jd1-jd0);
yr = iyr + (jd-jd0)./(jd1-jd0);

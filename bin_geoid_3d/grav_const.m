% GRAV_CONST: gravitational constants Re Rekm mu mukm gm J2=j2

% Ales Bezdek, bezdek@asu.cas.cz, 1/2015

Re=6378137;
Re=6378136.3;   %EGM2008, ECSS 2000
Rekm=6378.1363;   %EGM2008, ECSS 2000

mu=3.986004415e14;  %EGM2008, ECSS 2000
gm=3.986004415e14;  %EGM2008, ECSS 2000


j2=1082.6267e-6;   %EGM96

mu_sq=sqrt(mu);
J2=j2;
mukm=mu/1e9;

Ce=2*pi*Re; %equatorial circumference
Cekm=2*pi*Rekm;

mu_wgs84= 398600.5*1e9;            %// in m3 / s2

we=2*pi/86164.1;                  % angular velocity of the earth

function [cnm, snm, header, modelname]=cti_icgem_model(filename,nmax,t_yr,varargin)
% [cnm snm header modelname]=cti_icgem_model(filename,nmax,t_yr)
%     t_yr..epoch (years with decimals)
%     Nacte casove zavisle modely z icgem, museji byt predem pripraveny pomoci icgem2mat*
%     t_yr<0..nechci casovou zavislost
%     modely jsem prevedl na svuj format pomoci: icgem2mat*.m
% [cnm snm header modelname]=cti_icgem_model(filename,nmax,t_yr,'seasonal',0) 
%     bez sezonni slozky, jen stredni pole + trend

% 5/2018 pridavam clear nize, jinak si to zapamatovalo tvg cleny pro stat model
% 10/2017 predelano znacne: jediny load, varargin, snad to pocita dobre
%    uchovavam 
% 11/3/2015 neni to dobre, ted zbytecne mockrat nacitam model, chce to upravit
% AB, 5/3/12


global nmax2 filename2 t_yr2 seasonal2
global cnm2 snm2 header2 modelname2 n_t0 n_trnd n_acos n_asin cnm_t0 cnm_trnd snm_trnd cnm_acos snm_acos cnm_asin snm_asin yrs y1
% ecnm esnm: zatim odtud neexportuju
global cnm0 snm0

if isempty(nmax2); nmax2=-1; end
if isempty(t_yr2); t_yr2=-1; end

if nargin<2
   s=load(filename,'header');
   nmax=s.header.max_degree;
end

% Parsing optional arguments
seasonal=1;
i=1;
while i<=length(varargin)
   switch lower(varargin{i})
      case 'seasonal'
         seasonal=varargin{i+1};
   end
   i=i+2;
end


if ~strcmp(filename,filename2) || nmax2~=nmax
   n_t0=0;
   load(filename);
   filename2=filename;
   %    nmax2=nmax;
   %    cnm2=cnm; snm2=snm;
   header2=header;
   modelname2=modelname;
   nmax2=-1;
   fprintf('  loaded: %s\n',filename);
else
   cnm=cnm2;
   snm=snm2;
   header=header2;
   modelname=modelname2;
end

nmax1=nmax+1;
ii=1:nmax1;

if ~n_t0  %v modelu neni tvg
%    if nmax<header.max_degree
   if nmax2~=nmax
      if nmax<header.max_degree 
         cnm=cnm(ii,ii); snm=snm(ii,ii);
      end
      cnm2=cnm; snm2=snm;
      nmax2=nmax;
   end
   return
end

if  nargin<3
   t_yr=-1;
end

if t_yr<0   %nechci casovou zavislost
   if ~isempty(cnm0) && nmax2==nmax %mohu chtit najednou pouze stredni pole, kdyz uz jsem predtim chtel tvg
      cnm=cnm0; snm=snm0;
      cnm2=cnm; snm2=snm;
   elseif nmax2~=nmax
      if nmax<header.max_degree 
         cnm=cnm(ii,ii); snm=snm(ii,ii);
      end
      cnm2=cnm; snm2=snm;
      nmax2=nmax;
   end
   t_yr2=-1;
   return
end

%% v modelu je tvg
% nutno uschovat jeho puvodni shc (cnm0/snm0) pro tvg vypocet mezi volanimi procedury
if ~strcmp(filename,filename2) || nmax2~=nmax || t_yr2<0
   cnm=cnm(ii,ii,:); snm=snm(ii,ii,:);
   cnm0=cnm; snm0=snm;
   % zde orezu vsechny cleny na nmax2
   if n_t0  %zdali je v modelu tvg; da-li chybu, nebyl model zpracovan na tvg
      if n_trnd
         i1=cnm_trnd(:,1,1)<=nmax1 & cnm_trnd(:,2,1)<=nmax1;
%          n_trnd=sum(i1);
         cnm_trnd(~i1,:,:)=[];
         i1=snm_trnd(:,1,1)<=nmax1 & snm_trnd(:,2,1)<=nmax1;
         n_trnd=sum(i1);
         snm_trnd(~i1,:,:)=[];
      end
      if n_t0
         i1=cnm_t0(:,1,1)<=nmax1 & cnm_t0(:,2,1)<=nmax1;
         n_t0=sum(i1);
         cnm_t0(~i1,:,:)=[];
      end
      if n_acos
         i1=cnm_acos(:,1,1)<=nmax1 & cnm_acos(:,2,1)<=nmax1;
         n_acos=sum(i1);
         cnm_acos(~i1,:,:)=[];
         snm_acos(~i1,:,:)=[];
      end
      if n_asin
         i1=cnm_asin(:,1,1)<=nmax1 & cnm_asin(:,2,1)<=nmax1;
         n_asin=sum(i1);
         cnm_asin(~i1,:,:)=[];
         snm_asin(~i1,:,:)=[];
      end
   end
   nmax2=nmax;
   t_yr2=-1;
end

ys = [];
if t_yr2~=t_yr || seasonal2~=seasonal
   cnm=cnm0; snm=snm0;
   % seradit jednotlive roky a podle nich indexovat treti rozmer
   yrs = sort(unique(cnm_t0(:,3,:)));
   if yrs(1) == 0, yrs(1) = []; end  % prvni 'rok' muze byt nula, smazat
   y = 1;
   while t_yr > yrs(y)
    if y == length(yrs)
       y = y + 1;
       break
    end
    y = y + 1;
   end
   y1 = y-1;  % hledany index je o 1 mensi
   for i=1:n_trnd
      y = y1;
      % v pripade, ze nejsou vsechna data dostupna pro referencni epochu,
      % pouzijou se drivejsi
      while cnm_trnd(i,1,y) == 0
          y = y-1;
      end
      ys = [ys,y];
      n1=cnm_trnd(i,1,y);
      m1=cnm_trnd(i,2,y);
      i1= cnm_t0(:,1,y)==n1 & cnm_t0(:,2,y)==m1;
      t0_yr=cnm_t0(i1,3,y);
      trnd=cnm_trnd(i,3,y);
      dcnm=trnd*(t_yr-t0_yr);
      cnm_v(n1,m1)=cnm(n1,m1,y)+dcnm;
      trnd=snm_trnd(i,3,y);
      dsnm=trnd*(t_yr-t0_yr);
      snm_v(n1,m1)=snm(n1,m1,y)+dsnm;
   end
   if seasonal
      for i=1:n_acos
         y = y1;
         while cnm_acos(i,1,y) == 0
            y = y-1;
         end
         ys = [ys,y];
         n1=cnm_acos(i,1,y);
         m1=cnm_acos(i,2,y);
         i1= cnm_t0(:,1,y)==n1 & cnm_t0(:,2,y)==m1;
         t0_yr=cnm_t0(i1,3,y);

         cos1=cnm_acos(i,3,y);
         per=cnm_acos(i,4,y);
         dcnm=cos1*cos(2*pi/per*(t_yr-t0_yr));
         cnm_v(n1,m1)=cnm_v(n1,m1)+dcnm;

         cos1=snm_acos(i,3,y);
         per=snm_acos(i,4,y);
         dsnm=cos1*cos(2*pi/per*(t_yr-t0_yr));
         snm_v(n1,m1)=snm_v(n1,m1)+dsnm;

         sin1=cnm_asin(i,3,y);
         per=cnm_asin(i,4,y);
         dcnm=sin1*sin(2*pi/per*(t_yr-t0_yr));
         cnm_v(n1,m1)=cnm_v(n1,m1)+dcnm;

         sin1=snm_asin(i,3,y);
         per=snm_asin(i,4,y);
         dsnm=sin1*sin(2*pi/per*(t_yr-t0_yr));
         snm_v(n1,m1)=snm_v(n1,m1)+dsnm;
      end
   else
      fprintf('  no seasonal components\n');
   end
   t_yr2=t_yr;
   cnm2 = cnm; snm2 = snm;
   cnm = cnm_v; snm = snm_v;
   seasonal2=seasonal;
end
ys = unique(ys); % vsechny pouzite epochy
fprintf('  %s: time variability terms computed.\n',filename);
yr=floor(t_yr);
doy=365.25*(t_yr-yr);
dtn=doy2dtn(yr,doy+1); % pridam pul den, aby mi 1/1/2005 nedalo 31/12/2004
fprintf('     approximate epoch: %s\n',datestr(dtn));

% referencni epocha
yrr = floor(yrs(y1));
doyr = 365.25*(yrs(y1)-yrr);
dtnr = doy2dtn(yrr,doyr+1);
fprintf('     referential epoch: %s\n',datestr(dtnr));

function Pnm=alfs_Pnm_column(nmax,m,colatd)
% PNM=ALFS_PNM_COLUMN(NMAX,M,COLATD)  Associated Legendre functions of the first kind
%
% Standard forward column methods, HF, p. 281
% Input: nmax=maximum degree and order, m order, colatitude (rad)
% For given order m the function returns all the Pnm values with n=m:nmax.
% E.g. for m=0, nmax=N: P00 P10 P20 P30 etc PN0
%
% References:
% (HF) Holmes Featherstone 2002 A unified approach to the Clenshaw summation
%      and the recursive computation of very high degree and order normalised associated Legendre functions

% 20/7/12 Ales Bezdek, bezdek@asu.cas.cz
%% znovu
Ni=length(colatd);
Nj=nmax-m+1;
Pnm=zeros(Ni,Nj);

colat=colatd(:)/rad;
u=sin(colat);
t=cos(colat);

% sectoral seed, Eq 13
%Pnm(0,0)=0;
if m==0
   Pnm(:,1)=1;
elseif m==1
   Pnm(:,1)=sqrt(3)*u;
elseif m>1
   i=2*(2:m);
   i1=sqrt((i+ones(size(i)))./i);
   Pnm(:,1)=u.^m*sqrt(3)*prod(i1);
end

if m==nmax
   return
end

% recursion for non-sectoral Pnm: nmax==m+1
n=m+1;
anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
Pnm(:,2)=anm*t.*Pnm(:,1);

if m+1==nmax
   return
end

j=3;
for n=m+2:nmax
   anm=sqrt((2*n-1)*(2*n+1)/((n-m)*(n+m)));
   bnm=sqrt((2*n+1)*(n+m-1)*(n-m-1)/((n-m)*(n+m)*(2*n-3)));
   Pnm(:,j)=anm*t.*Pnm(:,j-1)-bnm*Pnm(:,j-2);
   j=j+1;
end

end

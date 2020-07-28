function cnm_U=cnm_normal(nmax)
%cnm_U=cnm_normal(nmax)
%normal field
%GRS80 suggests the following values:
%GM =3986005108 m3=s2; R = a = 6378137 m; C0;0 = 1;

% AB: 2/5/11
cnm_U=zeros(nmax+1,nmax+1);
cnm_U(1,1)=1;
cnm_U(3,1)=-4.841669e-4; 
cnm_U(5,1)=7.9030407e-7;
cnm_U(7,1)=1.687251e-9; 
cnm_U(9,1)=3.46098e-12;
cnm_U=cnm_U(1:nmax+1,1:nmax+1);

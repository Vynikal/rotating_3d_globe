% ICGEM2MAT   Reads geopotential coefficients from an ICGEM file and saves them in a mat file.
%
% Usage: 
% 
%       icgem2mat
% 
% finds all the ICGEM files (*.gfc) in the current directory,
% reads the geopotential coefficients, transforms them into Matlab variables:
%       header...structure with Icgem header information
%       cnm(n+1,m+1), snm(n+1,m+1)...harmonic coefficients C(n,m), S(n,m)
% 
% The new mat file with the same name is moved into 'data_icgem' subdirectory;
% the original gfc file is moved into 'data_icgem/gfc/' subdirectory.
% 
% Add the 'data_icgem' folder into your Matlab path.
% The model coefficients are then loaded by typing, e.g.:
% 
%          load egm2008
% 
% To display the C(2,0) zonal term type
% 
%          cnm(3,1)
% 
%
% See also compute_geopot_grids

% Ales Bezdek, bezdek@asu.cas.cz, 11/2012

clear
NMAX=360;
NMAX=1e100;  %it is possible to limit the maximum degree read from the gfc file 
adr_data='./';
adr_kam='./data_icgem/';

seznam_soub=dir(adr_data);
soub={seznam_soub.name};   %cell with filenames
for i=1:length(soub)
   jm=soub{i};
   if length(jm)>4 && strcmpi(jm(end-3:end),'.gfc')
      soub1=jm(1:end-4);
      fprintf('Gfc file processed: %s\n',soub1);
      filename=[adr_data soub1 '.gfc'];

      % Read header
      fid=fopen(filename);
      modelname=''; GM=0; ae=0; Lmax=0; errors=''; norm=''; tide='';

      s=fgets(fid);
      while(strncmp(s, 'end_of_head', 11) == 0 && sum(s)>=0)
         if (strncmp(s, 'product_type', 12)), product_type=strtrim(s(13:end)); end;
         if (strncmp(s, 'modelname', 9)), modelname=strtrim(s(10:end)); end;
         if (strncmp(s, 'earth_gravity_constant', 22)), GM=str2double(s(23:end)); end;
         if (strncmp(s, 'radius', 6)), ae=str2double(s(7:end)); end;
         if (strncmp(s, 'max_degree', 10)), Lmax=str2double(s(11:end)); end;
         if (strncmp(s, 'errors', 6)), errors=strtrim(s(7:end)); end;
         if (strncmp(s, 'norm', 4)), norm=strtrim(s(5:end)); end;
         if (strncmp(s, 'tide_system', 11)), tide=strtrim(s(12:end)); end;
         s=fgets(fid);
      end
      if sum(s)<0
         error_ab('Problem with reading the gfc file.')
      end

      header=struct('product_type',product_type,'modelname',modelname,'earth_gravity_constant',GM,'radius',ae,'max_degree',Lmax,'errors',errors,'norm',norm,'tide_system',tide);

      % read coefficients
      cnm=zeros(Lmax+1);
      snm=zeros(Lmax+1);
      ecnm=zeros(Lmax+1);
      esnm=zeros(Lmax+1);

      i_t0=0;
      i_trnd=0; %pocet clenu s trendem
      i_acos=0; %pocet clenu 
      i_asin=0; %pocet clenu 
      i_gfc=0;
      cnm_t0=[]; cnm_trnd=[]; snm_trnd=[]; cnm_acos=[]; snm_acos=[]; cnm_asin=[]; snm_asin=[]; 
      
      s=fgets(fid);
      while (length(s)>3)
         x=str2num(s(5:end));
         n=x(1)+1;
         m=x(2)+1;
         if n>NMAX || m>NMAX
            s=fgets(fid);
            continue;
         end
%         disp(s(1:4))
         if strcmp(s(1:4),'gfct')
            if isempty(cnm_t0)
               [status, result] =system(['grep -c gfct ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with t0'); end
               cnm_t0=zeros(i1,3);
            end
            i_t0=i_t0+1;
            cnm(n,m)=x(3);
            snm(n,m)=x(4);
%            if strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated')
               [yr,mn,dy]=ymd2cal(x(end)/1e4);
               yrd=jd2yr(cal2jd(yr,mn,dy));
               cnm_t0(i_t0,:)=[n m yrd];
%             elseif strcmp(header.errors,'calibrated_and_formal')
%             elseif strcmp(header.errors,'no')
%             end
            if (strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated') || strcmp(header.errors,'calibrated_and_formal')),
               ecnm(n,m)=x(5);
               esnm(n,m)=x(6);
            end
         elseif strcmp(s(1:3),'gfc')
            cnm(n,m)=x(3);
            snm(n,m)=x(4);
            if (strcmp(header.errors, 'formal') || strcmp(header.errors,'calibrated') || strcmp(header.errors,'calibrated_and_formal')),
               ecnm(n,m)=x(5);
               esnm(n,m)=x(6);
            end
            i_gfc=i_gfc+1;
         elseif strcmp(s(1:4),'trnd') || strcmp(s(1:3),'dot') 
            if isempty(cnm_trnd)
               [status, result] =system(['grep -c trnd ' filename]);
               i1=str2double(result);
               if i1==0; [status, result] =system(['grep -c dot ' filename]); i1=str2num(result); end
               if i1==0; error_ab('Problem with trnd'); end
               cnm_trnd=zeros(i1,3); snm_trnd=cnm_trnd;
            end
            i_trnd=i_trnd+1;
            cnm_trnd(i_trnd,:)=[n m x(3)];
            snm_trnd(i_trnd,:)=[n m x(4)];
         elseif strcmp(s(1:4),'acos')
            if isempty(cnm_acos)
               [status, result] =system(['grep -c acos ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with acos'); end
               cnm_acos=zeros(i1,4); snm_acos=cnm_acos;
            end
            i_acos=i_acos+1;
            cnm_acos(i_acos,:)=[n m x(3) x(end)];
            snm_acos(i_acos,:)=[n m x(4) x(end)];
         elseif strcmp(s(1:4),'asin')
            if isempty(cnm_asin)
               [status, result] =system(['grep -c asin ' filename]);
               i1=str2double(result);
               if i1==0; error_ab('Problem with asin'); end
               cnm_asin=zeros(i1,4); snm_asin=cnm_asin;
            end
            i_asin=i_asin+1;
            cnm_asin(i_asin,:)=[n m x(3) x(end)];
            snm_asin(i_asin,:)=[n m x(4) x(end)];
         else
            error_ab('A problem occured in gfc data.');
         end
         s=fgets(fid);
      end
      fclose(fid);

      modelname=header.modelname;

      %it is possible to limit the maximum degree read from the gfc file 
      n_gfc=i_gfc; n_t0=i_t0; n_trnd=i_trnd; n_acos=i_acos; n_asin=i_asin; 
      if n_t0; cnm_t0=cnm_t0(1:n_t0,:); end
      if n_trnd; cnm_trnd=cnm_trnd(1:n_trnd,:); snm_trnd=snm_trnd(1:n_trnd,:); end
      if n_acos; cnm_acos=cnm_acos(1:n_acos,:); snm_acos=snm_acos(1:n_acos,:); end
      if n_asin; cnm_asin=cnm_asin(1:n_asin,:); snm_asin=snm_asin(1:n_asin,:); end
      if n_t0~=n_trnd || n_acos~=n_asin
         error_ab('Problem with numbers of TVG terms.');
      end
%       fprintf('   gfc terms: %d, gfct: %d, trnd: %d, acos: %d, asin: %d\n',[n_gfc n_t0 n_trnd n_acos n_asin]);

      if ~exist(adr_kam,'file'), mkdir(adr_kam); end
      if ~exist([adr_kam 'gfc'],'file'), mkdir([adr_kam 'gfc']); end
      eval(sprintf('save %s%s.mat cnm snm ecnm esnm header modelname n_t0 n_trnd n_acos n_asin cnm_t0 cnm_trnd snm_trnd cnm_acos snm_acos cnm_asin snm_asin;',adr_kam,soub1));
      movefile([adr_data soub1 '.gfc'],[adr_kam 'gfc']);
      fprintf('  Resulting file %s.mat was moved into folder: %s\n',soub1,adr_kam);
      fprintf('  Original file  %s.gfc was moved into folder: %sgfc\n',soub1,adr_kam);
   end
end

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
      form=''; modelname=''; GM=0; ae=0; Lmax=0; errors=''; norm=''; tide='';

      s=fgets(fid);
      while(strncmp(s, 'end_of_head', 11) == 0 && sum(s)>=0)
         if (strncmp(s, 'format', 6)), form=strtrim(s(7:end)); end
         if (strncmp(s, 'product_type', 12)), product_type=strtrim(s(13:end)); end
         if (strncmp(s, 'modelname', 9)), modelname=strtrim(s(10:end)); end
         if (strncmp(s, 'earth_gravity_constant', 22)), GM=str2double(s(23:end)); end
         if (strncmp(s, 'radius', 6)), ae=str2double(s(7:end)); end
         if (strncmp(s, 'max_degree', 10)), Lmax=str2double(s(11:end)); end
         if (strncmp(s, 'errors', 6)), errors=strtrim(s(7:end)); end
         if (strncmp(s, 'norm', 4)), norm=strtrim(s(5:end)); end
         if (strncmp(s, 'tide_system', 11)), tide=strtrim(s(12:end)); end
         s=fgets(fid);
      end
      if sum(s)<0
         error_ab('Problem with reading the gfc file.')
      end

      header=struct('product_type',product_type,'modelname',modelname,'earth_gravity_constant',GM,'radius',ae,'max_degree',Lmax,'errors',errors,'norm',norm,'tide_system',tide,'format',form);

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
      cnm_t0=[]; cnm_trnd=[]; snm_trnd=[]; cnm_acos=[]; snm_acos=[]; cnm_asin=[]; snm_asin=[]; id = 1;
      line1 = ftell(fid); % zapamatovat posledni radku
      
    if form == 'icgem2.0' %jestli format je v2, projit vsechna data a vypsat datumy
        s=fgets(fid);
        while (length(s)>3)
            x=str2num(s(5:end));
            if ~strcmp(s(1:4),'gfc ')
                i_t0 = i_t0+1;
                y(i_t0) = x(7);
            end
            s=fgets(fid);
        end
        y = unique(y); y = sort(y); id = 1:length(y); T = containers.Map(y,id); % mapa s indexy/datumy
        i_t0=0; fseek(fid,line1,'bof'); % navrat na end_of_head
    end
    
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
         if strcmp(s(1:4),'gfc ')
            cnm(n,m,:)=x(3);
            snm(n,m,:)=x(4);
            if contains(header.errors,'calibrated') || contains(header.errors,'formal')
               ecnm(n,m)=x(5);
               esnm(n,m)=x(6);
            end
            i_gfc=i_gfc+1;
         elseif strcmp(s(1:4),'gfct')
            if isempty(cnm_t0)
               [status, result] = grep('-c','gfct  ', filename);
               i1=str2double(result.result);
               if i1==0; error_ab('Problem with t0'); end
               cnm_t0=zeros(i1,3);
            end
            if isempty(form) % pro stary format pouze jeden cas
                time = x(end);
                ti = 1;
            else
                time = x(end-1);
                ti = T(time); % vyhleda v mape index pro dany rok
            end
            i_t0=i_t0+1;
            cnm_t0(i_t0,:)=[n m time];
            cnm(n,m,ti)=x(3);
            snm(n,m,ti)=x(4);
            if contains(header.errors,'calibrated') || contains(header.errors,'formal')
               ecnm(n,m,ti)=x(5);
               esnm(n,m,ti)=x(6);
            end
         elseif strcmp(s(1:4),'trnd') || strcmp(s(1:3),'dot') 
            if isempty(cnm_trnd)
               [status, result] = grep('-c','trnd  ', filename);
               i1=str2double(result.result);
               if i1==0; [status, result] = grep('-c','dot  ', filename); i1=str2double(result.result); end
               if i1==0; error_ab('Problem with trnd'); end
               cnm_trnd=zeros(i1,3); snm_trnd=cnm_trnd;
            end
            i_trnd=i_trnd+1;
            if form == 'icgem2.0'
                t = x(end-1);
                t_trnd(i_trnd) = T(t);
            end
            cnm_trnd(i_trnd,:)=[n m x(3)];
            snm_trnd(i_trnd,:)=[n m x(4)];
         elseif strcmp(s(1:4),'acos')
            if isempty(cnm_acos)
               [status, result] = grep('-c','acos  ', filename);
               i1=str2double(result.result);
               if i1==0; error_ab('Problem with acos'); end
               cnm_acos=zeros(i1,4); snm_acos=cnm_acos;
            end
            i_acos=i_acos+1;
            if form == 'icgem2.0'
                t = x(end-2);
                t_acos(i_acos) = T(t);
            end
            cnm_acos(i_acos,:)=[n m x(3) x(end)];
            snm_acos(i_acos,:)=[n m x(4) x(end)];
         elseif strcmp(s(1:4),'asin')
            if isempty(cnm_asin)
               [status, result] = grep('-c','asin  ', filename);
               i1=str2double(result.result);
               if i1==0; error_ab('Problem with asin'); end
               cnm_asin=zeros(i1,4); snm_asin=cnm_asin;
            end
            i_asin=i_asin+1;
            if form == 'icgem2.0'
                t = x(end-2);
                t_asin(i_asin) = T(t);
            end
            cnm_asin(i_asin,:)=[n m x(3) x(end)];
            snm_asin(i_asin,:)=[n m x(4) x(end)];
         else
            error_ab('A problem occured in gfc data.');
         end
         s=fgets(fid);
      end
      fclose(fid);
      
      % pro format 2.0 preorganizovat matice do zadouciho tvaru
      % v t_trnd/t_acos... jsou indexy tretiho rozmeru zavisle na epose
      if form == 'icgem2.0'
        cnm_tr = zeros(length(unique(cnm_trnd(:,1)*100+cnm_trnd(:,2))),3,max(id)); snm_tr = cnm_tr; cnm_t = cnm_tr;
        cnm_ac = zeros(length(unique(cnm_acos(:,1)*100+cnm_acos(:,2)+cnm_acos(:,4)))-2,4,max(id)); cnm_as = cnm_ac; snm_ac = cnm_ac; snm_as = cnm_ac;
        j = 1;
        cnm_tr(j,:,t_trnd(1)) = cnm_trnd(1,:);
        snm_tr(j,:,t_trnd(1)) = snm_trnd(1,:);
        for i = 2:length(t_trnd)
            if cnm_trnd(i-1,1) ~= cnm_trnd(i,1) || cnm_trnd(i-1,2) ~= cnm_trnd(i,2)
                j = j+1;
            end
            cnm_tr(j,:,t_trnd(i)) = cnm_trnd(i,:);
            snm_tr(j,:,t_trnd(i)) = snm_trnd(i,:);
        end

        j = 1;
        cnm_t(j,:,t_trnd(1)) = cnm_t0(1,:);
        for i = 2:length(t_trnd)
            if cnm_t0(i-1,1) ~= cnm_t0(i,1) || cnm_t0(i-1,2) ~= cnm_t0(i,2)
                j = j+1;
            end
            cnm_t(j,:,t_trnd(i)) = cnm_t0(i,:);
        end
        j = 1;
        cnm_ac(j,:,t_acos(1)) = cnm_acos(1,:);
        cnm_ac(j+1,:,t_acos(1)) = cnm_acos(2,:);
        snm_ac(j,:,t_acos(1)) = snm_acos(1,:);
        snm_ac(j+1,:,t_acos(1)) = snm_acos(2,:);
        for i = 3:length(t_acos)
            if cnm_acos(i-1,1) ~= cnm_acos(i,1) || cnm_acos(i-1,2) ~= cnm_acos(i,2)
                j = j+2;
            end
            cnm_ac(j,:,t_acos(i)) = cnm_acos(i-1,:);
            cnm_ac(j+1,:,t_acos(i)) = cnm_acos(i,:);
            snm_ac(j,:,t_acos(i)) = snm_acos(i-1,:);
            snm_ac(j+1,:,t_acos(i)) = snm_acos(i,:);
        end
        j = 1;
        cnm_as(j,:,t_asin(1)) = cnm_asin(1,:);
        cnm_as(j+1,:,t_asin(1)) = cnm_asin(2,:);
        snm_as(j,:,t_asin(1)) = snm_asin(1,:);
        snm_as(j+1,:,t_asin(1)) = snm_asin(2,:);
        for i = 3:length(t_asin)
            if cnm_asin(i-1,1) ~= cnm_asin(i,1) || cnm_asin(i-1,2) ~= cnm_asin(i,2)
                j = j+2;
            end
            cnm_as(j,:,t_asin(i)) = cnm_asin(i-1,:);
            cnm_as(j+1,:,t_asin(i)) = cnm_asin(i,:);
            snm_as(j,:,t_asin(i)) = snm_asin(i-1,:);
            snm_as(j+1,:,t_asin(i)) = snm_asin(i,:);
        end
        cnm_trnd = cnm_tr; cnm_acos = cnm_ac; cnm_asin = cnm_as; snm_trnd = snm_tr; snm_acos = snm_ac; snm_asin = snm_as; cnm_t0 = cnm_t;
      end

      modelname=header.modelname;

      n_gfc=i_gfc; n_t0=i_t0; n_trnd=i_trnd; n_acos=i_acos; n_asin=i_asin; 
%       %it is possible to limit the maximum degree read from the gfc file 
%       if n_t0; cnm_t0=cnm_t0(1:n_t0,:,:); end
%       if n_trnd; cnm_trnd=cnm_trnd(1:n_trnd,:,:); snm_trnd=snm_trnd(1:n_trnd,:); end
%       if n_acos; cnm_acos=cnm_acos(1:n_acos,:,:); snm_acos=snm_acos(1:n_acos,:); end
%       if n_asin; cnm_asin=cnm_asin(1:n_asin,:,:); snm_asin=snm_asin(1:n_asin,:); end
      if n_t0~=n_trnd || n_acos~=n_asin
         error_ab('Problem with numbers of TVG terms.');
      end
%       fprintf('   gfc terms: %d, gfct: %d, trnd: %d, acos: %d, asin: %d\n',[n_gfc n_t0 n_trnd n_acos n_asin]);

      if ~exist(adr_kam,'file'), mkdir(adr_kam); end
      if ~exist([adr_kam 'gfc'],'file'), mkdir([adr_kam 'gfc']); end
      eval(sprintf('save %s%s.mat cnm snm ecnm esnm header modelname n_t0 n_trnd n_acos n_asin cnm_t0 cnm_trnd snm_trnd cnm_acos snm_acos cnm_asin snm_asin;',adr_kam,soub1));
      copyfile([adr_data soub1 '.gfc'],[adr_kam 'gfc']);
      fprintf('  Resulting file %s.mat was moved into folder: %s\n',soub1,adr_kam);
      fprintf('  Original file  %s.gfc was moved into folder: %sgfc\n',soub1,adr_kam);
   end
end

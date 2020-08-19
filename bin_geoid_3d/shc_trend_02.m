% m-file: shc_trend_02
% namalujeme si sekularni trend pro jeden rok
% nutno uzit cyklus for, protoze cti_icgem_model neni vektorova funkce

clear all

% ref_grav_model_name='eigen-6s';
ref_grav_model_name='EIGEN-6S2';
n = 5; m = 0;
nmax=5;
cnm1=zeros(24,nmax+1,nmax+1);
snm1=cnm1; cnm2=cnm1; snm2=cnm1;
dtn=zeros(24,1);
j=0;

for y=2000:2014
   for mo=1:12
       for d = 1:5:30
      j=j+1;
      
      epoch_yr=jd2yr(cal2jd(y,mo,d));
      dtn(j)=datenum(y,mo,d);
      ymd=y*10000+mo*100+d;
      [cnm, snm, header]=cti_icgem_model(lower(ref_grav_model_name),nmax,ymd);
      model=header.modelname;
      cnm1(j,:,:)=cnm;
      snm1(j,:,:)=snm;
      
      [cnm_noSeas, snm_noSeas]=cti_icgem_model(lower(ref_grav_model_name),nmax,ymd,'seasonal',0);
      % odecteme slozku konst+trend, zbude nam sinusova tvg, rijnova tvg, ta je treba na obr. 8 v HP
      cnm2(j,:,:)=cnm_noSeas;
      snm2(j,:,:)=snm_noSeas;
       end
   end
end
%% obrazek
clf
plot(dtn,cnm1(:,n+1,m+1)-0.9572e-6,'r.-')
hold on
plot(dtn,cnm2(:,n+1,m+1)-0.9572e-6,'.-')
title(sprintf('Závislost koeficientu C%d%d na èase podle modelu %s',n,m,ref_grav_model_name))
date = [];
for i = 1:length(dtn)
    if strncmp(datestr(dtn(i)),'01-Jan',6)
        date = [date;dtn(i)];
    end
end
xticks(date)
datetick('x','keepticks');
grid
legend 'konst+trend+sezonni' 'konst+trend'
% title(sprintf('%s: koeficient C(3,0)',lomitka(upper(ref_grav_model_name))));

grid on

w=1200; h=400; set(gcf, 'WindowStyle','normal','Position',[1,1,w,h],'PaperPositionMode','auto', 'PaperUnits','points','PaperSize',[w,h]);
soub = sprintf('shc_trend_02_%s_C%d%d.png',ref_grav_model_name,n,m);
eval(sprintf('print -dpng %s',soub))

crop(soub)

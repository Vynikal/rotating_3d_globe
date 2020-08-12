% m-file: shc_trend_02
% namalujeme si sekularni trend pro jeden rok
% nutno uzit cyklus for, protoze cti_icgem_model neni vektorova funkce

clear all

% ref_grav_model_name='eigen-6s';
ref_grav_model_name='EIGEN-6S4_v2';

nmax=3;
cnm1=zeros(24,nmax+1,nmax+1);
snm1=cnm1; cnm2=cnm1; snm2=cnm1;
dtn=zeros(24,1);
j=0;

for y=2000:2014
   for m=1:12
       for d = 1:5:30
      j=j+1;
      
      epoch_yr=jd2yr(cal2jd(y,m,d));
      dtn(j)=datenum(y,m,d);
      ymd=y*10000+m*100+d;
      [cnm, snm, header]=cti_icgem_model(lower(ref_grav_model_name),nmax,ymd,'period',0.5);
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
plot(dtn,cnm1(:,4,1)-0.9572e-6,'r.-')
hold on
plot(dtn,cnm2(:,4,1)-0.9572e-6,'.-')
datetick('x');
grid 
legend 'konst+trend+sezonni' 'konst+trend'
% title(sprintf('%s: koeficient C(3,0)',lomitka(upper(ref_grav_model_name))));
set(gca,'ytick',2000:2015)
grid on

w=1200; h=400; set(gcf, 'WindowStyle','normal','Position',[1,1,w,h],'PaperPositionMode','auto', 'PaperUnits','points','PaperSize',[w,h]);
eval(sprintf('print -dpng -r0 shc_trend_02_%s.png',ref_grav_model_name));

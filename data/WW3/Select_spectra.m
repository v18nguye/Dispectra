clear all
close all
clc

str_date = '201001';
pnt_name = 'NODE008919';

% spectrum file
SPC_file = sprintf('ww3.%s_%s_spec.nc', pnt_name, str_date);
% image saved file name
sfn = sprintf('wind_sea_%s_%s.png',pnt_name,str_date);

% corresponding integrated wave parameters timeseries file
IWP_file = sprintf('ww3.%s_%s_para.nc', pnt_name, str_date);

% load files
SPC = Read_NetCDF(SPC_file);
IWP = Read_NetCDF(IWP_file);

% mask for found wave system(s)
mask_p0 = ~isnan(IWP.phs0); % partition 0 is for wind sea
mask_p1 = ~isnan(IWP.phs1); % partition 1 is for most energetic swell
mask_p2 = ~isnan(IWP.phs2); % partition 2 is for second most energetic swell
mask_p3 = ~isnan(IWP.phs3); % ....

%
% plot first spectrum with only 1 wave system found (WS)
%

% select date with detected wind sea, without swell
MatTime_WS = IWP.MatTime(mask_p0 & ~mask_p1 & ~mask_p2 & ~mask_p3);

% select date of first spectrum
MatTime = MatTime_WS(1);

% time indices in files
b1 = SPC.MatTime;
b2 = IWP.MatTime;
i1 = find(abs(SPC.MatTime-MatTime) < 1e-10);
i2 = find(abs(IWP.MatTime-MatTime) < 1e-10);

% get spectrum
d = SPC.direction([1:end,1]);
freq  = SPC.frequency;
ang = SPC.direction([1:end,1]);
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% polar plot of the spectrum
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);
pcolor(fx,fy,Efth)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

% add detail on detected partitions

annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '1 system found :', ...
    sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2))})

saveas(gcf,sfn);
sfn = sprintf('wind_sea2_%s_%s.png',pnt_name,str_date);

% Create a second graph for segmentation by using the Watersed method
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)));
fx2(2:26,2:31) = fx;
fy2(2:26,2:31) = fy; 
Efth2(2:26,2:31) = Efth;
pcolor(fx2,fy2,Efth2);
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '1 system found :', ...
    sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2))})
saveas(gcf,sfn)
sfn = sprintf('1_swell_%s_%s.png',pnt_name,str_date);

%
% plot first spectrum with only 1 wave system found (S1)
%

% select date with detected swell (wind sea), without wind sea (swell)
MatTime_S1 = IWP.MatTime(~mask_p0 & mask_p1 & ~mask_p2 & ~mask_p3);

% select date of first spectrum
MatTime = MatTime_S1(1);

% time indices in files
i1 = find(abs(SPC.MatTime-MatTime) < 1e-10);
i2 = find(abs(IWP.MatTime-MatTime) < 1e-10);

% get spectrum
freq  = SPC.frequency;
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% polar plot of the spectrum
tit = sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime));
figure('Name',tit)
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);
pcolor(fx,fy,Efth)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)

% add detail on detected partitions
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '1 system found :', ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2))})

saveas(gcf,sfn)
sfn = sprintf('1_swell2_%s_%s.png',pnt_name,str_date);


% Create a second graph for segmentation by using the Watersed method
tit = sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime));
figure('Name',tit);
fx2(2:26,2:31) = fx;
fy2(2:26,2:31) = fy; 
Efth2(2:26,2:31) = Efth;
pcolor(fx2,fy2,Efth2);
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '1 system found :', ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2))})
saveas(gcf,sfn)
sfn = sprintf('wind_sea_swell_%s_%s.png',pnt_name,str_date);


%
% plot first spectrum with 2 wave systems (SW + S1) found
%

% select date with detected wind sea + 1 swell
MatTime_WS_S1 = IWP.MatTime(mask_p0 & mask_p1 & ~mask_p2 & ~mask_p3);

% select date of first spectrum
MatTime = MatTime_WS_S1(1);

% time indices in files
i1 = find(abs(SPC.MatTime-MatTime) < 1e-10);
i2 = find(abs(IWP.MatTime-MatTime) < 1e-10);

% get spectrum
freq  = SPC.frequency;
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% polar plot of the spectrum
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);
pcolor(fx,fy,Efth)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)

% add detail on detected partitions
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '2 systems found :', ...
    sprintf('Hs(SW) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2))})
saveas(gcf,sfn)
sfn = sprintf('wind_sea_swell2_%s_%s.png',pnt_name,str_date);

% Create a second graph for segmentation by using the Watersed method
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)));
fx2(2:26,2:31) = fx;
fy2(2:26,2:31) = fy; 
Efth2(2:26,2:31) = Efth;
pcolor(fx2,fy2,Efth2);
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '2 systems found :', ...
    sprintf('Hs(SW) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2))})
saveas(gcf,sfn)
sfn = sprintf('wind_sea_2swells_%s_%s.png',pnt_name,str_date);

%
% plot first spectrum with 3 wave systems found
%

% select date with detected wind sea + 2 (1) swell
MatTime_WS_S1_S2 = IWP.MatTime(mask_p0 & mask_p1 & mask_p2 & ~mask_p3);

% select date of first spectrum
MatTime = MatTime_WS_S1_S2(1); 

% time indices in files
i1 = find(abs(SPC.MatTime-MatTime) < 1e-10);
i2 = find(abs(IWP.MatTime-MatTime) < 1e-10);

% get spectrum
freq  = SPC.frequency;
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% polar plot of the spectrum
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);
pcolor(fx,fy,Efth)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)

% add detail on detected partitions
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '3 systems found :', ...
    sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2)), ...
    sprintf('Hs(S2) = %4.1f m ; Dir(S2) =%3d deg',IWP.phs2(i2),IWP.pdir2(i2))})

saveas(gcf,sfn)
sfn = sprintf('wind_sea_2swells2_%s_%s.png',pnt_name,str_date);

% Create a second graph for segmentation by using the Watersed method
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)));
fx2(2:26,2:31) = fx;
fy2(2:26,2:31) = fy; 
Efth2(2:26,2:31) = Efth;
pcolor(fx2,fy2,Efth2);
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(tit,'fontsize',10)
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '3 systems found :', ...
    sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2)), ...
    sprintf('Hs(S2) = %4.1f m ; Dir(S2) =%3d deg',IWP.phs2(i2),IWP.pdir2(i2))})
saveas(gcf,sfn)
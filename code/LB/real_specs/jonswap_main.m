clear all
close all
clc

rng(1)
%%
%load needed paths.
addpath(genpath('./simu'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/data/WW3'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/blasso'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/lasso'))

addpath(genpath('E:/IMT/intern2019/data/WW3'))
addpath(genpath('E:/IMT/intern2019/code/LB/bla  sso'))
addpath(genpath('E:/IMT/intern2019/code/LB/lasso'))
%%
% spec file parameters
str_date = '201001';
pnt_name = 'NODE008919';

% spectrum file
SPC_file = sprintf('ww3.%s_%s_spec.nc', pnt_name, str_date);

% corresponding integrated wave parameters timeseries file
IWP_file = sprintf('ww3.%s_%s_para.nc', pnt_name, str_date);

%%
% load files
SPC = Read_NetCDF(SPC_file);
IWP = Read_NetCDF(IWP_file);

% mask for found wave system(s)
mask_p0 = ~isnan(IWP.phs0); % partition 0 is for wind sea
mask_p1 = ~isnan(IWP.phs1); % partition 1 is for most energetic swell
mask_p2 = ~isnan(IWP.phs2); % partition 2 is for second most energetic swell
mask_p3 = ~isnan(IWP.phs3); % ....

%% select the first spectrum with only one wave system (WS)

% select date with detected wind sea, without swell
MatTime_WS = IWP.MatTime(mask_p0 & mask_p1 & mask_p2 & ~mask_p3);

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
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
%theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% convert to hte polar coordinate.
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);

%%
% SFW

% range of each parameter.
    %   range(:,:,1) =[Hmin Tmin ; Hmax Tmax]
    %   range(:,:,2) =[cmin theta0min; cmax theta0max]
range(:,:,1) = [0.1 10; 2 30];
range(:,:,2) = [10 0.01*pi; 30 2*pi];

% the JONSWAP shape's parameter.
gam = 3.3;

% number of parameter elements.
    % N = [N_H, N_T, N_c, N_theta0]
N = [15 15 15 18];

y = reshape(Efth,[],1); % spec observation.

% simulate.
s = jonswap_simu('jonswap', gam, range, freq, theta);

% swf method simulation parameters.
opts.param_grid = s.test_grid(N); % create a parameter grid
opts.A = s.atom(opts.param_grid);
opts.atom = s.atom;
opts.datom = s.datom;
opts.B = s.range;
opts.cplx = s.cplx;
lambda_lambdaMax = .01;
lambdaMax = norm(opts.A'*y,inf);
opts.lambda = lambda_lambdaMax*lambdaMax;
opts.maxIter = 100;
opts.tol = 1.e-5;
opts.disp = true;
opts.mergeStep = 0.05; %0.01

% resolve the system.
tic
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW4d( y , opts );
toc

y1_reconstruct = opts.atom(param_SFW_blasso)*x_SFW_blasso;

%%

%
% polar plot of the original occean wave spectrum
%
figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
pcolor(fx,fy,Efth)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

% add detail on detected partitions
annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '3 systems found :', ...
    sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2)), ...
    sprintf('Hs(S2) = %4.1f m ; Dir(S2) =%3d deg',IWP.phs2(i2),IWP.pdir2(i2))})


%
% polar plot of four spec having the largest coefficients and a spec by suming
% all of them.
%

x_SFW_blasso_abs = abs(x_SFW_blasso);
v = zeros(length(x_SFW_blasso),4);
i = zeros(1,4);
% find the largest coefficient.
[~,i1] = max(x_SFW_blasso_abs);
m1 = x_SFW_blasso(i1,1);
% assign the i1th value to the minimum value of the vector
x_SFW_blasso_abs(i1) = min(x_SFW_blasso_abs);
v(i1,1) = m1;
i(1,1) = i1;

% find the second largest coefficient.
[~,i2] = max(x_SFW_blasso_abs);
m2 = x_SFW_blasso(i2,1);
%...
x_SFW_blasso_abs(i2) = min(x_SFW_blasso_abs);
v(i2,2) = m2;
i(1,2) = i2;
% find the third largest coefficient.
[~,i3] = max(x_SFW_blasso_abs);
m3 = x_SFW_blasso(i3,1);
%...
x_SFW_blasso_abs(i3) = min(x_SFW_blasso_abs);
v(i3,3) = m3;
i(1,3) = i3;
% find the fourth largest coefficient.
[~,i4] = max(x_SFW_blasso_abs);
m4 = x_SFW_blasso(i4,1);
%...
x_SFW_blasso_abs(i4) = min(x_SFW_blasso_abs);
v(i4,4) = m4;
i(1,4) = i4;

% plot separately four spec elements
figure('Name',sprintf('Spectrum Elements'))
for k =1:4
subplot(2,2,k)
pcolor(fx,fy,reshape(abs(opts.atom(param_SFW_blasso)*v(:,k)),size(Efth)))
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title(['H_T[',num2str(param_SFW_blasso(1,i(1,k))),' ',num2str(param_SFW_blasso(2,i(1,k))),'] c_theta[',num2str(param_SFW_blasso(3,i(1,k))),' ',num2str(param_SFW_blasso(4,i(1,k))),'] - a: ', num2str(x_SFW_blasso(i(1,k),1))]); 
end

% plot the sum of four specs
figure('Name',sprintf('Spectrum recovered by the SFW method'))
pcolor(s.fx,s.fy,reshape(real(y1_reconstruct),size(Efth)))
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title('By all spectrum units')

figure('Name',sprintf('Spectrum recovered by the SFW method'))
y2_reconstruct =0;
for k = 1:4
    y2_reconstruct = y2_reconstruct + opts.atom(param_SFW_blasso)*v(:,k);
end
pcolor(s.fx,s.fy,reshape(real(y2_reconstruct),size(Efth)))
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
title('By all four largest coef spec units')
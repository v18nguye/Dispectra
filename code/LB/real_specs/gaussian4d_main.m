clear all
close all
clc

%%
%load needed paths.
addpath(genpath('./simu'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/data/WW3'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/blasso'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/lasso'))

addpath(genpath('E:/IMT/intern2019/data/WW3'))
addpath(genpath('E:/IMT/intern2019/code/LB/blasso'))
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
MatTime_WS = IWP.MatTime(mask_p0 & mask_p1 & ~mask_p2 & ~mask_p3);

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
Efth = SPC.efth([1:end,1],:,i1);

% convert to hte polar coordinate.
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);

%% sfw method

% initiate a simulation for the continuous dict case (blasso).
rangexy = [2*min(min(fx)) 2*min(min(fy)); 2*max(max(fx)) 2*max(max(fy))];
range_sigxy2 = [0.001 0.001;2 2];
b_range = cat(3,rangexy,range_sigxy2);
N = 30;
y = reshape(Efth,[],1); % spec observation.
simu_opts = gaussian_4d_simu('4dgaussian', b_range, fx, fy);

% parameters
opts.param_grid = simu_opts.test_grid(N);
opts.A = simu_opts.atom(opts.param_grid);

opts.atom = simu_opts.atom;
opts.datom = simu_opts.datom;
opts.B = simu_opts.p_range;
opts.cplx = simu_opts.cplx;

lambda_lambdaMax = .2;
lambdaMax = norm(opts.A'*y,inf);

opts.lambda = lambda_lambdaMax*lambdaMax;
opts.maxIter = 200;
opts.tol = 1.e-5;
opts.disp = true;
tic
opts.mergeStep = .1;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW4d( y , opts );
toc
y_blasso = opts.atom(param_SFW_blasso)*x_SFW_blasso;

%%
if true 
    % polar plot of the spectrum
    figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
    pcolor(fx,fy,Efth)
    shading flat
    cb = colorbar;
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

    % add detail on detected partitions
    annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
    'backgroundcolor','w',...
    'string',{...
    '2 systems found :', ...
    sprintf('Hs(SW) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2)), ...
    sprintf('Hs(S1) = %4.1f m ; Dir(S1) =%3d deg',IWP.phs1(i2),IWP.pdir1(i2))}) 
    
    figure('Name',sprintf('WS recovered by SFW for %s (%s)', pnt_name, datestr(MatTime)))
    pcolor(fx,fy,reshape(y_blasso,size(Efth)))
    shading flat
    cb = colorbar;
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
end

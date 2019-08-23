clear all
close all
clc

%%
%load needed paths.
addpath(genpath('./simu'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/data/WW3'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/blasso'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/lasso'))

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

%%
% SFW

% range of each parameter.
    %   range(:,:,1) =[Hmin Tmin ; Hmax Tmax]
    %   range(:,:,2) =[cmin theta0min; cmax theta0max]
range(:,:,1) = [0.1 1; 15 20];
range(:,:,2) = [1 0; 20 2*pi];

% the JONSWAP shape's parameter.
gam = 3.3;

% number of parameter elements.
    % N = [N_H, N_T, N_c, N_theta0]
N = [20 20 20 36];

y = reshape(Efth,[],1); % spec observation.

% simulate.
simu_opts = jonswap_simu('jonswap', gam, range);

% swf method simulation parameters.
opts.param_grid = simu_opts.test_grid(N); % create a parameter grid
opts.A = simu_opts.atom(opts.param_grid);
opts.atom = simu_opts.atom;
opts.datom = simu_opts.datom;
opts.B = simu_opts.range;
opts.cplx = simu_opts.cplx;
lambda_lambdaMax = .01;
lambdaMax = norm(opts.A'*y,inf);
opts.lambda = lambda_lambdaMax*lambdaMax;
opts.maxIter = 200;
opts.tol = 1.e-5;
opts.disp = true;
opts.mergeStep = .01;

% resolve the system.
tic
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW4d( y , opts );
toc

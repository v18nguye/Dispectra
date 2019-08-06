clear all
close all
clc

%%
%load needed paths.
addpath(genpath('./simu'))
addpath(genpath('./lasso'))
addpath(genpath('./dictgen'))
addpath(genpath('./blasso'))

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
theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);

% convert to hte polar coordinate.
[ffreq,ttheta] = meshgrid(freq,theta);
[fx,fy] = pol2cart(ttheta,ffreq);

%% sfw method

% initiate a simulation in the continuous dict case (blasso).
b_range = [1.5*min(min(fx)) 1.5*min(min(fy)); 1.5*max(max(fx)) 1.5*max(max(fy))];
N = 100;
y = reshape(Efth,[],1);
simu_opts = real_spec_blasso('2dgaussian', b_range, fx, fy);

% parameters
optsb.param_grid = simu_opts.test_grid(N);
optsb.A = simu_opts.atom(optsb.param_grid);

optsb.atom = simu_opts.atom;
optsb.datom = simu_opts.datom;
optsb.B = simu_opts.p_range;
optsb.cplx = simu_opts.cplx;

lambda_lambdaMax = .01;
lambdaMax = norm(optsb.A'*y,inf);

optsb.lambda = lambda_lambdaMax*lambdaMax;
optsb.maxIter = 10000;
optsb.tol = 1.e-5;
optsb.disp = true;

optsb.mergeStep = .01;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW( y , optsb );

%%
if true 
    % polar plot of the spectrum
    figure('Name',sprintf('Wave Spectrum for %s (%s)', pnt_name, datestr(MatTime)))
    pcolor(fx,fy,Efth)
    shading interp
    cb = colorbar;
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

    % add detail on detected partitions
    annotation('textbox',[0.05 0.20 0.01 0.01],'FitBoxToText','on',...
        'backgroundcolor','w',...
        'string',{...
        '1 system found :', ...
        sprintf('Hs(WS) = %4.1f m ; Dir(WS) =%3d deg',IWP.phs0(i2),IWP.pdir0(i2))})
end


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
% initiate parameters.
K = 3;      % number of sources in the simulated signal.
SNR = inf;  % input snr.
rng(1)      % set the seed.
N = 100; % number of sampling coordination on each axis.
rangeuxy = [-4 -4; 4 4]; % range of uxy
rangesigxy = [0.1 0.1; 1 1];
range = cat(3,rangeuxy,rangesigxy);
r_spec = 6; % the radius of spec
s_spec = 0.3; % the sampling of spec

%%
% radius
radi = 0:s_spec:r_spec;

% angle
theta = 0:pi/18:2*pi;

[rradi, ttheta] = meshgrid(radi, theta);

%
% convert the polar coordinate elements to the cartesian coordinate.
%

[fx,fy] = pol2cart(ttheta,rradi);

simu_opts = real_spec_blasso('2dgaussian',range,fx,fy);    % select the type of data to be simulated. possible choices are: 'doa' , 'gaussian' , 'dgaussian')

[ paramr , coef , y] = simu_opts.simu(K,SNR);  % simulate the observation vector

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
optsb.maxIter = 100;
optsb.tol = 1.e-5;
optsb.disp = true;

optsb.mergeStep = .01;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = real_spec_SFW( y , optsb );
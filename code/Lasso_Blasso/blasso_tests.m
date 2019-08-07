clear
close all;
clc

addpath(genpath('./simu'))
addpath(genpath('./lasso'))
addpath(genpath('./dictgen'))
addpath(genpath('./blasso'))

N = 100;    % size of the dictionary A
K = 5;      % number of sources in the simulated signal
SNR = inf;  % input snr
rng(1)      % set the seed
range = [-4 -4; 4 4]; % range of uxy
r_spec = 6; % the radius of spec
s_spec = 0.1; % the sampling of spec

%%
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%
simu_opts = blasso_simu('2dgaussian',range,r_spec,s_spec);    % select the type of data to be simulated. possible choices are: 'doa' , 'gaussian' , 'dgaussian')

[ paramr , coef , y] = simu_opts.simu(K,SNR);  % simulate the observation vector


%%
%%%%%%%%%%%%%%%%%%%%%%
% Methods parameters %
%%%%%%%%%%%%%%%%%%%%%%
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

%%
%%%%%%%%%%%%%%
% SFW-blasso %
%%%%%%%%%%%%%%
optsb.mergeStep = .01;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW( y , optsb );


%%
%%%%%%%%%%%%%%
% Figures    %
%%%%%%%%%%%%%%
figure
scatter3(paramr(1,:),paramr(2,:),coef)
xlim(simu_opts.p_range(:,1))
ylim(simu_opts.p_range(:,2))
zlim([(min(coef)-2) (max(coef)+2)])
hold all
scatter3(param_SFW_blasso(1,:),param_SFW_blasso(2,:),x_SFW_blasso)


figure
plot(fc_SFW_lasso)
xlabel('iter.')
ylabel('lasso cost-function')


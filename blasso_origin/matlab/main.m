clear
close all
clc
addpath(genpath('./fcts'))
addpath(genpath('./simu'))

M = 101;    % size of the observation vector y : M*1
N = 100;    % size of the dictionary A : M*N
K = 2;      % number of sources in the simulated signal
SNR = inf;  % input snr
rng(1)      % set the seed

%%
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%
simu_opts = simu_prop( 'gaussian' , M );    % select the type of data to be simulated. possible choices are: 'doa' , 'gaussian' , 'dgaussian')

[ theta_y , x_y , y ] = simu_opts.simu(K,SNR);  % simulate the observation vector


%%
%%%%%%%%%%%%%%%%%%%%%%
% Methods parameters %
%%%%%%%%%%%%%%%%%%%%%%
opts.param_grid = simu_opts.test_grid(N);
opts.A = simu_opts.atom(opts.param_grid);
A = opts.A;
opts.atom = simu_opts.atom;
opts.datom = simu_opts.datom;
opts.B = simu_opts.p_range;
opts.cplx = simu_opts.cplx;
%opts.m = simu_opts.m;

lambda_lambdaMax = .01;
lambdaMax = norm(opts.A'*y,inf);

opts.lambda = lambda_lambdaMax*lambdaMax;
opts.maxIter = 1000;
opts.tol = 1.e-5;
opts.disp = true;

%%
%%%%%%%%%
% Fista %
%%%%%%%%%
opts.L = max(eig(opts.A'*opts.A));
opts.xinit = zeros(N,1);
[x_FISTA_lasso , fc_FISTA_lasso , fc_FISTA_lassodual] = fista( y , opts );


%%
%%%%%%%%%%%%%
% FW-blasso %
%%%%%%%%%%%%%
[param_FW_blasso, x_FW_blasso , fc_FW_blasso , fc_FW_lasso , fc_FW_lassodual ] = FW( y , opts );


%%
%%%%%%%%%%%%%%
% SFW-blasso %
%%%%%%%%%%%%%%
opts.mergeStep = .1;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual, A_sfw ] = SFW( y , opts );



%%
%%%%%%%%%%%%%%
% Figures %
%%%%%%%%%%%%%%
figure
stem(theta_y,x_y)
xlim(simu_opts.p_range)
hold all
plot(opts.param_grid,x_FISTA_lasso)
plot(param_FW_blasso,x_FW_blasso,'xg')
plot(param_SFW_blasso,x_SFW_blasso,'xr')
legend('GT','Fista','F.-W.','Sliding-F.-W.')

figure
semilogy(fc_FISTA_lasso)
hold all
plot(fc_FW_lasso)
plot(fc_SFW_lasso)
xlabel('iter.')
ylabel('lasso cost-function')
legend('Fista','F.-W.','Sliding-F.-W.')
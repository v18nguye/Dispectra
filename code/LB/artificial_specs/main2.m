%
% This code simulate the SFW method used to segment the artificial spectra
% modeled by the gaussian atom where their four parameters vary.
%
clear all
close all
clc

%%
% load needed paths.

addpath(genpath('./continu_dict/simu'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/blasso'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/lasso'))

%%
% initiate the spec's parameters.

K = 3;      % number of sources in the simulated signal.
SNR = inf;  % input snr.
rng(1)      % set the seed.
N = 10; % number of sampling coordination on each axis.
rangeuxy = [-4 -4; 4 4]; % range of uxy
rangesigxy = [0.1 0.1; 1 1];
range = cat(3,rangeuxy,rangesigxy);
r_spec = 6; % the radius of spec
s_spec = 0.3; % the sampling of spec

% radius of the spectra
radi = 0:s_spec:r_spec;

% range of angle in the polar coordinate
theta = 0:pi/18:2*pi;

[rradi, ttheta] = meshgrid(radi, theta);

% convert the polar coordinate elements to the cartesian coordinate.
[fx,fy] = pol2cart(ttheta,rradi);

%%
% simulation the spec.

simu_opts = gaussian_4d_simu('4dgaussian',range,fx,fy);    % select the type of data to be simulated. possible choices are: 'doa' , 'gaussian' , 'dgaussian')

[ paramgt , coef , y] = simu_opts.simu(K,SNR);  % simulate the observation vector

%%
% reconstruct the spec by using the SFW method. 

% method's parameters
opts.param_grid = simu_opts.test_grid(N);
opts.A = simu_opts.atom(opts.param_grid);

opts.atom = simu_opts.atom;
opts.datom = simu_opts.datom;
opts.B = simu_opts.p_range;
opts.cplx = simu_opts.cplx;

lambda_lambdaMax = .1;
lambdaMax = norm(opts.A'*y,inf);

opts.lambda = lambda_lambdaMax*lambdaMax;
opts.maxIter = 200;
opts.tol = 1.e-5;
opts.disp = true;

tic
opts.mergeStep = .1;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW4d( y , opts );
toc

% the reconstructed spectrum
y_rec = opts.atom(param_SFW_blasso)*x_SFW_blasso;

%%
% visualization

% draw the spectrum entering to the algorithm.
if true
    figure('Name','the simulated spectrum')
    s= pcolor(fx,fy,reshape(y,size(fx)));
    shading flat
    
    cb = colorbar();
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

end

% draw the spectrum units which form the spectrum.
if false
    figure('Name','constructing spectrum units')
    for i = 1:length(coef)
        subplot(ceil(length(coef)/2),2,i);
        pcolor(fx,fy,reshape(simu_opts.atom(paramgt(:,i))*coef(i,1),size(fx)))
        shading flat
        title(['u_xy[',num2str(paramgt(1,i)),' ',num2str(paramgt(2,i)),'] sig2_xy[',num2str(paramgt(3,i)),' ',num2str(paramgt(4,i)),'] - a: ', num2str(coef(i,1))]); 

    end
end

% draw wave spectrum units obtained by SWF method.
if false
    figure('Name','spectrum units reconstructed by SFW')
    for i = 1: length(x_SFW_blasso)
        subplot(ceil(length(x_SFW_blasso)/2),2,i)
        pcolor(fx,fy,reshape(simu_opts.atom(param_SFW_blasso(:,i))*x_SFW_blasso(i,1),size(fx)))
        shading flat
        title(['u_xy[',num2str(param_SFW_blasso(1,i)),' ',num2str(param_SFW_blasso(2,i)),'] sig2_xy[',num2str(param_SFW_blasso(3,i)),' ',num2str(param_SFW_blasso(4,i)),'] - a: ', num2str(x_SFW_blasso(i,1))]); 
    end
end

% draw the reconstructed spectrum obtained by SFW method.
if true
    figure('Name','the reconstructed spectrum')
    s= pcolor(fx,fy,reshape(y_rec,size(fx)));
    shading flat
    
    cb = colorbar();
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

end
clear all
close all
clc

%%
%load needed paths.
addpath(genpath('./continu_dict/simu'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/blasso'))
addpath(genpath('/homes/v18nguye/Documents/intern2019/code/LB/lasso'))

%%
% initiate parameters.
K = 2;      % number of sources in the simulated signal.
SNR = inf;  % input snr.
%rng(1)      % set the seed.
N = 10; % number of sampling coordination on each axis.
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
%%
% simulation
simu_opts = gaussian_4d_simu('4dgaussian',range,fx,fy);    % select the type of data to be simulated. possible choices are: 'doa' , 'gaussian' , 'dgaussian')

[ paramgt , coef , y] = simu_opts.simu(K,SNR);  % simulate the observation vector

% draw the simulated spec.
if true
    
    %figure('Name','simulated spectra')
    figure
    s= pcolor(fx,fy,reshape(y,size(fx)));
    %title('3 separated spectra units')
    shading interp
    
    cb = colorbar('Direction', 'normal');
    set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')
    colormap parula

end

% draw the spectrum units
if true
    figure
    for i = 1:length(coef)
        subplot(ceil(length(coef)/2),2,i);
        pcolor(fx,fy,reshape( simu_opts.atom(paramgt(:,i))*coef(i,1),size(fx)))
        shading interp
        title(['u_xy[',num2str(paramgt(1,i)),' ',num2str(paramgt(2,i)),'] sig2_xy[',num2str(paramgt(3,i)),' ',num2str(paramgt(4,i)),'] - a: ', num2str(coef(i,1))]); 

    end
end

%%
% Use the SFW to solve the blasso problem
if true
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
    optsb.maxIter = 200;
    optsb.tol = 1.e-5;
    optsb.disp = true;

    tic
    optsb.mergeStep = .01;
    [param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW4d( y , optsb );
    toc

    x_SFW_blasso2 = [];
    param_SFW_blasso2 = [];
end

% draw wave spectrum units obtained by SWF method.
if true
    
    figure('Name','SFW-Blasso');
    x_sfw_size = size(x_SFW_blasso);
    x_sfw_size = x_sfw_size(1);
    
    % Drawing the spectra units that have the coefficients large enough.
  
    for k = 1: x_sfw_size
        if abs(x_SFW_blasso(k,1)) > 0.1
            x_SFW_blasso2 = [x_SFW_blasso2; x_SFW_blasso(k,1)];
            param_SFW_blasso2 = [param_SFW_blasso2 param_SFW_blasso(:,k)];
        end
    end
    
    x_sfw_size2 = size(x_SFW_blasso2);
    x_sfw_size2 = x_sfw_size2(1);
    
    for k = 1: x_sfw_size2        
        subplot(ceil(x_sfw_size2/3),3,k);
        pcolor(fx,fy,reshape(optsb.atom(param_SFW_blasso2(:,k)),size(fx)));
        shading interp
        title(['u_xy[',num2str(param_SFW_blasso2(1,k)),' ',num2str(param_SFW_blasso2(2,k)),'] sig2_xy[',num2str(param_SFW_blasso2(3,k)),' ',num2str(param_SFW_blasso2(4,k)),'] - a: ', num2str(x_SFW_blasso2(k,1))]); 
    end
end
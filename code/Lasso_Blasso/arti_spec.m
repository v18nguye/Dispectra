clear
close all;
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

% simulated parameters.
dict_simu = false;
lasso_spec_simu = true;

% parameters in the discret dictionary case (lasso).
r_spec = 6; % the spec radius.
s_spec = 0.3; % the sampling step for the spec radius.
r = 4; % the radius of gaussian mean coordinations. 
s = 0.1; % the sampling step

% parameters in the continuous dictionary case (blasso).
N = 100; % number of sampling coordination on each axis.
b_range = []; % the gaussian mean coordination region.
if lasso_spec_simu
    % for guaranting that mean coordination blasso region covers the one of
    % lasso.
	b_range = [-r -r; r r]; 
else
    % for guaranting that the mean coordination lasso region covers the one
    % of blasso.
    b_range = [-r*(sqrt(2)/2) -r*(sqrt(2)/2); r*(sqrt(2)/2) r*(sqrt(2)/2)];
end
%%
% creat/load a discrete dictionary (lasso).

if dict_simu
% if true, creat dict.
    disp('Create a dictionary !');
    [dict, dez, xe, ye, u_cordi] =  dic_gen(r,s,r_spec,s_spec);
    % a structure saving variables.
    strt.dict = dict;
    strt.dez = dez;
    strt.xe = xe;
    strt.ye = ye;
    strt.u_cordi = u_cordi;
    save('dict_infos.mat','-struct', 'strt');
else
% if false, just load a existed dict.
    disp('Load a dictionary !');
    load('dict_infos.mat','dict','dez','xe','ye','u_cordi');
end

%%
% initiate a simulation in the continuous dict case (blasso).
simu_opts = blasso_simu('2dgaussian', b_range, r_spec, s_spec);

%%
% simulate a tested spec.
if lasso_spec_simu
    % if true, simulate a spec by using the discrete dict in the lasso case.
    disp('Simulate spec by using the discrete dict in the lasso case !');
    [y, pst, coef] = simu_spec(K, dict, SNR);
    paramgt = u_cordi(:,pst);
    
else
    % if false, simulate a spec by atom in the blasso case.
    disp('Simulate spec by using the atom in the blasso case !')
    [ paramgt , coef , y] = simu_opts.simu(K,SNR);
    
end

% draw the simulated spec.
if true
    
    figure
    pcolor(xe,ye,reshape(y,dez))
    title('simulated spectrum')
    shading interp
end

%%
% method parameters in the blasso case.
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

% method parmeters in the lasso case.
optsl.A = dict;
lambda_lambdaMax = .01;

lambdaMax = norm(optsl.A'*y,inf);

optsl.lambda = lambda_lambdaMax*lambdaMax;
optsl.maxIter = 10000;
optsl.tol = 1.e-5;
optsl.disp = true;

%%
%%%%%%%%%
% Fista %
%%%%%%%%%
dict_siz = size(dict);
N = dict_siz(2); % number of the dictionary elements. 
optsl.L = max(eig(optsl.A'*optsl.A));
optsl.xinit = zeros(N,1);

optsl.L = max(eig(optsl.A'*optsl.A));
tic
[x_FISTA_lasso , fc_FISTA_lasso , fc_FISTA_lassodual] = fista( y , optsl );
toc

% draw spec units obtained by fista method.
if true
    
    % seuillage de x_FISTA_lasso
    for i = 1:N
        if abs(x_FISTA_lasso(i,1)) < 10^-1
            x_FISTA_lasso(i,1) = 0;
        end
    end
    % find the nonzero indices and their values
    [row_x,col_x,v_x] = find(x_FISTA_lasso);
    % the number of nonzero elements
    nonzero_eles = size(row_x);
    nonzero_eles = nonzero_eles(1);
    figure('Name','fista-lasso');
    for k =1:1:nonzero_eles

        % take each component indice and its value.
        fista_compo = zeros(N,1);
        fista_compo(row_x(k,1)) = v_x(k,1);
        fista_y = dict*fista_compo;

        subplot(ceil(nonzero_eles/3),3,k);
        pcolor(xe,ye,reshape(fista_y,dez));
        shading interp
        % uxy = [ux uy] - the mean values along 2 axis.
        % a - amplitude.
        title(['u_xy[',num2str(u_cordi(1,row_x(k,1))),'',num2str(u_cordi(2,row_x(k,1))),']', ' - a: ', num2str(v_x(k,1))]);    

    end
end

%%
%%%%%%%%%%%%%%
% SFW-blasso %
%%%%%%%%%%%%%%
tic
optsb.mergeStep = .01;
[param_SFW_blasso, x_SFW_blasso , fc_SFW_blasso , fc_SFW_lasso , fc_SFW_lassodual ] = SFW( y , optsb );
toc

% draw wave spec units obtained by SWF method.
if true
    
    figure('Name','SFW-Blasso');
    x_sfw_size = size(x_SFW_blasso);
    x_sfw_size = x_sfw_size(1);
    
    for k = 1: x_sfw_size
        
        subplot(ceil(x_sfw_size/3),3,k);
        pcolor(xe,ye,reshape(optsb.atom(param_SFW_blasso(:,k)),dez));
        shading interp
        % uxy = [ux uy] - the mean values along 2 axis.
        % a - amplitude.
        title(['u_xy[',num2str(param_SFW_blasso(1,k)),' ',num2str(param_SFW_blasso(2,k)),']', ' - a: ', num2str(x_SFW_blasso(k,1))]); 
    end
end

%%
%%%%%%%%%%%%%%
% Figures    %
%%%%%%%%%%%%%%

% draw ground truth, fista results, sfw results on 3D coordination.
% axis x representing the gaussian mean value along x axis.
% axis y representing the gaussian mean value along y axis.
% axis z representing the altitude of the gaussian wave specs.

if true    
    figure
    
    % draw the groud truth
    scatter3(paramgt(1,:),paramgt(2,:),coef)
    xlim(simu_opts.p_range(:,1))
    ylim(simu_opts.p_range(:,2))
    zlim([(min(coef)-2) (max(coef)+2)])
    hold all
    
    % draw the fista results
        % seuillage de x_FISTA_lasso
        for i = 1:N
            if abs(x_FISTA_lasso(i,1)) < 10^-1
                x_FISTA_lasso(i,1) = 0;
            end
        end
        % find the nonzero indices and their values
        [row_x,col_x,v_x] = find(x_FISTA_lasso);
        
    scatter3(param_SFW_blasso(1,:),param_SFW_blasso(2,:),v_x')
    
    % draw the sfw results
    scatter3(u_cordi(1,row_x(:,1)'),u_cordi(2,row_x(:,1)),x_SFW_blasso)

    
    figure
    plot(fc_SFW_lasso)
    xlabel('iter.')
    ylabel('lasso cost-function')
end





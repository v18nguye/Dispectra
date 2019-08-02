clear
close all;
clc

addpath(genpath('./simu'))
addpath(genpath('./lasso'))
addpath(genpath('./dictgen'))

snr = inf; % signal noise ratio input.
rng = 1; % set the seed.
N = 2; % number of used spectrum units
load('dict_infos.mat','dict','dez','xe','ye') % load the dictionary informations.

%%
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%

[y, pst, ampt] = simu_spec(N, dict, snr);

%%
%%%%%%%%%%%%%%%%%%%%%%
% Methods parameters %
%%%%%%%%%%%%%%%%%%%%%%

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
N = dict_siz(2); % number of the dictioinary elements 
optsl.L = max(eig(optsl.A'*optsl.A));
optsl.xinit = zeros(N,1);

[x_FISTA_lasso , fc_FISTA_lasso , fc_FISTA_lassodual] = fista( y , optsl );

% find the nonzero indices and their values
[row_x,col_x,v_x] = find(x_FISTA_lasso);

% the number of nonzero elements
nonzero_eles = size(row_x);
nonzero_eles = nonzero_eles(1);

%%
%%%%%%%%%%%%%%
% Figures    %
%%%%%%%%%%%%%%

figure
pcolor(xe,ye,reshape(y,dez))
title('simulated spectrum density')
shading interp
saveas(gcf, 'n5_simulated_spec.png');

% component size of each element in the dictionary
component_size = size(dict);
component_size = component_size(2);

ori_pst_sort = sort(pst);
ori_nonzeros_eles = size(ori_pst_sort);
ori_nonzeros_eles = ori_nonzeros_eles(2);

row_plot_max = max(nonzero_eles,ori_nonzeros_eles);

figure
% draw the ground truth spectrum density components.
for k = 1:1:ori_nonzeros_eles
    
    subplot(2,row_plot_max,k);
    pcolor(xe,ye,reshape(dict(:,ori_pst_sort(1,k)),dez));
    shading interp
    % o_i - orginal indice.
    % a - amplitude.
    title(['o_i: ', num2str(ori_pst_sort(1,k)), ' - a: ', num2str(ampt(k,1))]);    
end

% draw the spectrum density components found by FISTA.
for k =1:1:nonzero_eles
    
    % take each component indice and its value.
    fista_compo = zeros(component_size,1);
    fista_compo(row_x(k,1)) = v_x(k,1);
    fista_y = dict*fista_compo;
    
    subplot(2,row_plot_max,row_plot_max+k);
    pcolor(xe,ye,reshape(fista_y,dez));
    shading interp
    % f_i - fista indice of each component.
    % a - amplitude.
    title(['f_i: ', num2str(row_x(k,1)), ' - a: ', num2str(v_x(k,1))]);    
    
end
saveas(gcf, 'n5_spec_units.png');


figure
semilogy(fc_FISTA_lasso)
xlabel('iter.')
ylabel('lasso cost-function')
saveas(gcf, 'n5_cost_function.png');


figure
stem(pst,ampt,'LineStyle', 'none');
hold all
stem(row_x',v_x','LineStyle', 'none');
xlabel('indice');
ylabel('amplitude');
legend('Ground Truth','FISTA');
saveas(gcf, 'n5_indice_amplitude.png');



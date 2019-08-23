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
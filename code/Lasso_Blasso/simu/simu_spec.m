% simu_spec.m 
%
% This function generates a spectrum density y which are combined by several
% spectrum units being in the generated dictionary.
% 
% Usage: [ y, x_gt, r, ampt]  = simu_spec( k, A, SNR)
% 
% k - the number of pricipal components that form a signal.
%
% A - the dictionary.
%
% SNR - the signal noise ratio.
%
% Return: 
%
% y - the simulated spectrum.
%
% x_gt - the ground truth in the equation "y = Ax". 
%
% pst - component positions.
%
% ampt - amplitude of each corresponding component.
%
% Written by: Van-Khoa NGUYEN
% Email: van-khoa.nguyen@imt-atlantique.net
% Created: July 2019


function [ y, pst, ampt] = simu_spec( k, A, SNR)

%
% Here we simulate k spectrum units that are closed together. The maximum
% number of these spectrum units are 3.
%

% get the number of spectrum units.

matrix_siz = size(A);

num_of_spec = matrix_siz(2);

% generate a random number in the range of the number of spectrum units.

if (0 < k) && (k <= 3)
    
    r = randperm(num_of_spec,1);
    
    % vector contains at most the closed 3 spectrum unit positions. 
    r_vec = r:10:r+10*(k-1);
    r_vec = mod(r_vec,num_of_spec +1);
    
    % the ground truth vector.
    x = zeros(num_of_spec,1);
    
    x(r_vec) = randn(k,1);
    
    % generate the signal y.
    y0 = A*x;
    
    % add noise.
    n = randn(size(y0));
    
    n = n/norm(n)*norm(y0)*10^(-SNR/20); 
    
    y0 = y0 + n;
    
    y = y0;
    
    
    pst = r_vec;
    
    [~,~,val_x] = find(x);
    
    ampt = val_x;
    
     
else
    
    r = randperm(num_of_spec,2);
    
    % vector contains at most the closed 3 spectrum unit positions. 
    r_vec1 = r(1):1:r(1)+1*2;
    
    r_vec2 = r(2):1:r(2)+1*(k-4);
    
    r_vec = [r_vec1 r_vec2];
    r_vec = mod(r_vec,num_of_spec +1);
    
    r_vec = sort(r_vec);
    
    % the ground truth vector.
    x = zeros(num_of_spec,1);
    
    x(r_vec) = randn(k,1);
    
    % generate the signal y.
    y0 = A*x;
    
    % add noise
    n = randn(size(y0));
    
    n = n/norm(n)*norm(y0)*10^(-SNR/20);
    
    y0 = y0 + n;
    
    y= y0;
    
    pst = r_vec;
    
    [~,~,val_x] = find(x);
    
    ampt = val_x;
         
end

end
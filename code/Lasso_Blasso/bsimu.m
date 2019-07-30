function [ opts ] = bsimu( dico , m)
if(nargin==0)
    dico='2dgaussian';
    m=51;
elseif(nargin==1)
    m=51;
end

switch dico
    case '2dgaussian'
        p_range = [0;1]; % as r
        sigmax = 0.4;
        sigmay = 0.2;
        atom = @(param) atom_gaussian(param,m,sigmax,sigmay);
        datom = @(param) datom_gaussian(param,m,sigmax,sigmay);
        cplx = false;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
    case 'd2dgaussian'
        p_range = [0;1];
        sigmax = 0.4;
        sigmay = 0.2;
        atom = @(param) atom_dgaussian(param,m,sigmax); 
        datom = @(param) datom_dgaussian(param,m,sigmay);
        cplx = false;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
        
end
opts.cplx=cplx;
opts.p_range = p_range;
opts.atom = atom;
opts.datom = datom;
opts.simu =simu;
opts.m = m;
end


function [param,coeff,y] = signal(k,SNR,p_range,atom,cplx)
param = rand(1,k)*(p_range(2)-p_range(1))+p_range(1);
if(cplx)
    coeff = randn(k,1)+1i*randn(k,1);
else
    coeff = randn(k,1);
end
y = atom(param)*coeff;
n = randn(size(y));
n = n/norm(n)*norm(y)*10^(-SNR/20);
y=y+n;
end
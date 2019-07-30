function [ opts ] = simu_prop( dico , m)
if(nargin==0)
    dico='gaussian';
    m=51;
elseif(nargin==1)
    m=51;
end

switch dico
    case 'gaussian'
        p_range = [0;1]; % must to change this range for 2D Blasso.
        sigma = 5.e-2;
        atom = @(param) atom_gaussian(param,m,sigma);
        datom = @(param) datom_gaussian(param,m,sigma);
        cplx = true;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
    case 'dgaussian'
        p_range = [0;1];
        sigma = 5.e-2;
        atom = @(param) atom_dgaussian(param,m,sigma); 
        datom = @(param) datom_dgaussian(param,m,sigma);
        cplx = true;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
    case 'bigaussian'
        p_range = [0;1];
        sigma = 5.e-2;
        delta = .075;
        atom = @(param) atom_bigaussian(param,m,sigma,delta); 
        datom = @(param) datom_bigaussian(param,m,sigma,delta);
        cplx = false;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
    case 'doa'
        p_range = [ -pi/2 , pi/2 ];
        D = 1;
        W = 4;
        atom = @(param) atom_doa(param,m,D,W);
        datom = @(param) datom_doa(param,m,D,W);
        cplx = true;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
end
opts.cplx=cplx;
opts.p_range = p_range;
opts.atom = atom;
opts.datom = datom;
opts.simu =simu;
opts.test_grid = @(N) make_grid(N,p_range);
opts.m = m;
end



function A = atom_gaussian( param , m , sigma )
% m is the size of the observation vector.
t = linspace(0,1,m)';
T = repmat(t,1,length(param))-repmat(param,m,1);
A = exp(-0.5*T.^2/sigma^2);
A = A*diag(1./sqrt(sum(A.^2,1)));
end

function A = atom_dgaussian( param , m , sigma )
t = linspace(0,1,m)';
T = repmat(t,1,length(param))-repmat(param,m,1);
A = (T./sigma^2).*exp(-0.5*T.^2/sigma^2);
A = A*diag(1./sqrt(sum(A.^2,1)));
end



function A = atom_bigaussian( param , m , sigma , delta )
t = linspace(0,1,m)';
T1 = repmat(t,1,length(param))-repmat(param+delta,m,1);
T2 = repmat(t,1,length(param))-repmat(param-delta,m,1);
A = exp(-0.5*T1.^2/sigma^2)+exp(-0.5*T2.^2/sigma^2);
A = A*diag(1./sqrt(sum(A.^2,1)));
end


function A = datom_gaussian( param , m , sigma )
A = zeros(m,length(param));
for i=1:length(param)
    t = linspace(0,1,m)'-param(i);
    e = exp(-t.^2/(2*sigma^2));
    N = norm(e);
    A(:,i) = ((N*t/sigma^2).*e-((t'*e.^2)/(sigma^2*N))*e)/N^2;
end
end


function A = datom_dgaussian( param , m , sigma )
A = zeros(m,length(param));
for i=1:length(param)
    t = linspace(0,1,m)'-param(i);
    e = exp(-t.^2/(2*sigma^2));
    b = (t/sigma^2).*e;
    db = ((t/sigma^2).^2-1/sigma^2).*e;
    N = norm(b);
    A(:,i) = (db*N-(db'*b/N)*b)/N^2;
end
end



function A = datom_bigaussian( param , m , sigma , delta )
% % % TODO ! !
A = zeros(m,length(param));
for i=1:length(param)
    t = linspace(0,1,m)'-param(i);
    e = exp(-t.^2/(2*sigma^2));
    b = (t/sigma^2).*e;
    db = ((t/sigma^2).^2-1/sigma^2).*e;
    N = norm(b);
    A(:,i) = (db*N-(db'*b/N)*b)/N^2;
end
end

function A = atom_doa( param , m , D , W )
t = (-(m-1)/2:(m-1)/2)';
A =  1/sqrt(m)*exp(-2*1i*pi*D/W*t*sin(param));
end


function A = datom_doa( param , m , D , W )
t = (-(m-1)/2:(m-1)/2)';
A =  (-2*1i*pi*D/W*t*cos(param)).*(1/sqrt(m)*exp(-2*1i*pi*D/W*t*sin(param)));
end



function G = make_grid( N , B )
if(nargin==0)
    N=100;
end
G = B(1):(B(2)-B(1))/(N-1):B(2);
end


function [param,coeff,y, A_simu] = signal(k,SNR,p_range,atom,cplx)
param = (rand(1,k))*(p_range(2)-p_range(1))+p_range(1);
if(cplx)
    coeff = randn(k,1)+1i*randn(k,1);
else
    coeff = randn(k,1);
end

A_simu = atom(param);
y = atom(param)*coeff;
n = randn(size(y));
n = n/norm(n)*norm(y)*10^(-SNR/20);
y=y+n;
end
function [ opts ] = bsimu( dico)

switch dico
    case '2dgaussian'
        p_range = [0;5]; %  range of ux
        sigmax2 = 0.4;
        sigmay2 = 0.2;
        uy = 1;
        r = 7; % the radius of arti spec
        atom = @(ux) atom_2dgaussian(ux,r,sigmax2,sigmay2,uy);
        datom = @(ux) datom_2dgaussian(ux,r,sigmax2,sigmay2,uy);
        cplx = false;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
        
end
opts.cplx=cplx;
opts.p_range = p_range;
opts.atom = atom;
opts.datom = datom;
opts.simu =simu;
opts.test_grid = @(N) make_grid(N,p_range);
end



function A  = atom_2dgaussian(ux, r, sigx, sigy, uy)

n_ux = size(ux);
spec = [];

for uxi = 1:n_ux(2)
    % the mean values along 2 directions x,y
    m = [ux(1,uxi), uy]';

    % the deviation values along 2 directions x,y
    sig = [sigx, sigy]';

    % radius
    radi = 0:0.1:r;

    % angle
    theta = 0:pi/36:2*pi;

    [rradi, ttheta] = meshgrid(radi, theta);

    %
    % convert the polar coordinate elements to the cartesian coordinate.
    %

    [x,y] = pol2cart(ttheta,rradi);

    xy = cat(3,x,y);

    %
    % calculate the spectral density.
    %

    % subtract the mean values
    xym = bsxfun(@minus,xy,shiftdim(m,-2));

    % the reciprocal of sig times the xym
    rsigxym = bsxfun(@times,shiftdim(1./(sig),-2),xym);

    XY = sum(rsigxym.*rsigxym,3);

    spec0 = (1/(2*pi*nthroot(sigx*sigy,4)))*exp(-0.5*XY);

    spec0 = reshape(spec0,[],1);

    spec = [spec, spec0];

end

A = spec;
    
end


function A = datom_2dgaussian(ux, r, sigx, sigy, uy)

n_ux = size(ux);
dspec = [];

for uxi = 1:n_ux(2)
    % the mean values along 2 directions x,y
    m = [ux(1,uxi), uy]';

    % the deviation values along 2 directions x,y
    sig = [sigx, sigy]';

    % radius
    radi = 0:0.1:r;

    % angle
    theta = 0:pi/36:2*pi;

    [rradi, ttheta] = meshgrid(radi, theta);

    %
    % convert the polar coordinate elements to the cartesian coordinate.
    %

    [x,y] = pol2cart(ttheta,rradi);

    xy = cat(3,x,y);

    %
    % calculate the spectral density.
    %

    % subtract the mean values
    xym = bsxfun(@minus,xy,shiftdim(m,-2));

    % the reciprocal of sig times the xym
    rsigxym = bsxfun(@times,shiftdim(1./(sig),-2),xym);

    XY = sum(rsigxym.*rsigxym,3);

    spec0 = (1/(2*pi*nthroot(sigx*sigy,4)))*exp(-0.5*XY);
    
    duxi = (x -ux(1,uxi))*(1/(sigx));
    
    dspec0 = duxi.*spec0;

    dspec0 = reshape(dspec0,[],1);

    dspec = [dspec, dspec0];

end

A = dspec;

end

function G = make_grid( N , B )
if(nargin==0)
    N=100;
end
G = B(1):(B(2)-B(1))/(N-1):B(2);
end

function [param,coeff,y, A_simu] = signal(k,SNR,p_range,atom,cplx)
param = rand(1,k)*(p_range(2)-p_range(1))+p_range(1);
if(cplx)
    coeff = randn(k,1)+1i*randn(k,1);
else
    coeff = randn(k,1);
end
A_simu = atom(param );
y = atom(param)*coeff;
n = randn(size(y));
n = n/norm(n)*norm(y)*10^(-SNR/20);
y=y+n;
end
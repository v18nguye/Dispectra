function [ opts ] = blasso_simu( dico, range, r_spec,s_spec)

% disco - the type of atom.
% range - the 2D  gaussian mean coordination region.
% r - the radius of spec.


switch dico
    case '2dgaussian' % atom gaussian in two dimensions
        
        p_range = range; %  range of mean values in both axis (x,y)
        sigmax = 0.4;
        sigmay = 0.2;
        atom = @(uxy) atom_2dgaussian(uxy,r_spec,s_spec,sigmax,sigmay);
        datom = @(uxy) datom_2dgaussian(uxy,r_spec,s_spec,sigmax,sigmay);
        cplx = false;
        simu = @(k,SNR) signal(k,SNR,p_range,atom,cplx);
        
end
opts.cplx= cplx;
opts.p_range = p_range;
opts.atom = atom;
opts.datom = datom;
opts.simu = simu;
opts.test_grid = @(N) make_grid(N,p_range);
end


function A  = atom_2dgaussian(uxy, r,s, sigx, sigy)

n_uxy = size(uxy);
spec = [];

for ui = 1:n_uxy(2)
    
    % the mean values along 2 directions x,y
    ux = uxy(1,ui);
    uy = uxy(2,ui);
    
    m = [ux, uy]';

    % the deviation values along 2 directions x,y
    sig = [sigx, sigy]';

    % radius
    radi = 0:s:r;

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

    spec0 = (1/(2*pi*nthroot(sigx*sigy,2)))*exp(-0.5*XY);

    spec0 = reshape(spec0,[],1);

    spec = [spec, spec0];

end

A = spec;
    
end


function A = datom_2dgaussian(uxy, r,s, sigx, sigy)

n_uxy = size(uxy);
dspec = [];

for ui = 1:n_uxy(2)
    
    % the mean values along 2 directions x,y
    ux = uxy(1,ui);
    uy = uxy(2,ui);
    
    m = [ux, uy]';

    % the deviation values along 2 directions x,y
    sig = [sigx, sigy]';

    % radius
    radi = 0:s:r;

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

    spec0 = (1/(2*pi*nthroot(sigx*sigy,2)))*exp(-0.5*XY);
    
    dux = (x -ux)*(1/(sigx));
    duy = (y -uy)*(1/(sigy));
    
    dspec0x = reshape(dux.*spec0,[],1);
    
    dspec0y = reshape(duy.*spec0,[],1);

    dspec0 = [dspec0x dspec0y];

    dspec = [dspec, dspec0];

end

A = dspec;

end



function G = make_grid( N , B )

if(nargin==0)
    N=100;
end

Gx = B(1,1):(B(2,1)-B(1,1))/(N-1):B(2,1);
Gy = B(1,2):(B(2,2)-B(1,2))/(N-1):B(2,2);

[X,Y] = meshgrid(Gx, Gy);

G =[reshape(X,1,[]);reshape(Y,1,[])];

end


function [param,coeff,y] = signal(k,SNR,p_range,atom,cplx)

%paramx = rand(1,k)*(p_range(2,1)-p_range(1,1))+p_range(1,1);
%paramy = rand(1,k)*(p_range(2,2)-p_range(1,2))+p_range(1,2);
paramx = rand(1,k)*0.2*(p_range(2,1)-p_range(1,1))+p_range(1,1);
paramy = rand(1,k)*0.2*(p_range(2,2)-p_range(1,2))+p_range(1,2);
param = [paramx; paramy];

if(cplx)
    coeff = randn(k,1)+1i*randn(k,1);
else
    coeff = randn(k,1)+2;
end
y = atom(param)*coeff;
n = randn(size(y));
n = n/norm(n)*norm(y)*10^(-SNR/20);
y=y+n;
end
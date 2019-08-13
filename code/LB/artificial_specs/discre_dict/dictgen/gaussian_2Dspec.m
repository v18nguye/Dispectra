% gaussian_2Dspec.m
%
% This function aims to generate artificial spectrums which take a form of
% the gaussian distribution. The spectrums are drawn on the polar
% coordinate system.
%
% Usage:  spec = gaussian_spec( r, s, mx, my, sigx, sigy)
%
% s - sampling step
%
% r - the radius of a spectrum in the polar coordinate.
%
% mx - the mean of the 2D spec along the axis x.
%
% my - the mean of 2D spec along the axis y.
%
% sigx - the standard deviation along the axis x.
%   Default = 0.4.
%
% sigy - the standard deviation along the axis y.
%   Default = 0.2.
%
% Return:
%
% spec - a spec vector which represents the spectrum density.
%
% spec_size - the original size of the spec under the matrix form.
%
% xe - coordinations of dictionary element along to the x axis.
% 
% ye - coordinations of dictionary element along to the y axis.
%
% Written by: Van-Khoa NGUYEN
% Email: van-khoa.nguyen@imt-atlantique.net
% Created: July 2019

function [spec, spec_size, xe, ye] = gaussian_2Dspec( r, s, mx, my, sigx, sigy)

if (nargin < 5), sigx = 0.4; end
if (nargin < 6), sigy = 0.2; end

% the mean values along 2 directions x,y
m = [mx, my]';

% the deviation values along 2 directions x,y
sig = [sigx, sigy]';

%
% create the polar coordination.
%

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
rsigxym = bsxfun(@times,shiftdim(1./sig,-2),xym);

XY = sum(rsigxym.*rsigxym,3);

spec0 = (1/(2*pi*sqrt(sigx*sigy)))*exp(-0.5*XY);

%
% represent the spectrum on the cartesian coordinate.
%

%figure
%pcolor(x,y,spec0);
%shading interp

%
% reshape a spec matrix to a vector coresponding to an element of
% dictionary.
%

spec_size = size(spec0);

spec = reshape(spec0,[],1);

xe = x;

ye = y;

% reshape to the original shape by using 'reshape(spec,spec_size)'

end






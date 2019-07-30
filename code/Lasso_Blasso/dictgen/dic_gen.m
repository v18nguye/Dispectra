% dic_gen.m
%
% This function generates a dictionary which contains artificial spectral
% wave elements.
%
% Usage:  dict = dic_gen( r, s, re, se)
%
% r - the radius of dictionary which is used to put a limit to dictionary
% element coordinations.
%
% s - the sampling step
%
% re - the radius of dictionary elements which is used to generate
% spectrum units.
%
% se - the sampling step of dictionary elements
%
% Return:
%
% dict - a dictionary taking a matrix form whose elements are columns which represent pectrum units generated
% by gaussian_2D.
%
% dez - dictionary element size which is the original size of each spectrum
% element taking a matrix form.
%
% xe - coordinations of dictionary element along to the x axis.
% 
% ye - coordinations of dictionary element along to the y axis.
%
% Written by: Van-Khoa NGUYEN
% Email: van-khoa.nguyen@imt-atlantique.net
% Created: July 2019.

function [dict, dez, xe, ye] = dic_gen( r, s, re, se)

%
% generate the polar coordinations.
%

radi = 0:s:r;

thetha = 0:pi/18:2*pi;

[rradi, ttheta] = meshgrid(radi, thetha);

%
% convert to the polar coordinate elements to the cartesian coordinate.
%

[x,y] = pol2cart(ttheta, rradi);

%
% generate dictionary.
%

dict0 = [];

dez0 = [];

xe0 = [];

ye0 = [];

[kr, kc] = size(x);

for i = 1:kr
    
    for j = 1:kc
        
        % e represents for 'element'
        [dict_e, e_size, xxe, yye] = gaussian_2Dspec(re,se,x(i,j),y(i,j));
        
        dict0 = [dict0 dict_e];
        
        dez0 = e_size;
        
        xe0 = xxe;
        
        ye0 = yye;
    end
    
end

dict = dict0; 

dez = dez0;

xe = xe0;

ye = ye0;

end
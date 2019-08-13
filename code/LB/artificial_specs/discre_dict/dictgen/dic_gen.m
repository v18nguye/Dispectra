% dic_gen.m
%
% This function generates a dictionary which contains artificial spectral
% wave elements.
%
% Usage:  dict = dic_gen( r, s, r_spec, s_spec)
%
% r - the radius of mean coordinations.
%
% s - the sampling step for mean coordinations.
%
% r_spec - the radius of specs.
%
% s_spec - the sampling step for the spec radius.
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
% u_cordi - the mean coordinations corresponding to each dict element.
%
% Written by: Van-Khoa NGUYEN
% Email: van-khoa.nguyen@imt-atlantique.net
% Created: July 2019.

function [dict, dez, xe, ye, u_cordi] = dic_gen( r, s, r_spec, s_spec)

%
% generate the polar coordinations.
%

radi = 0:s:r;

thetha = 0:pi/36:2*pi;

[rradi, ttheta] = meshgrid(radi, thetha);

%
% convert to the polar coordinate elements to the cartesian coordinate.
%

[x,y] = pol2cart(ttheta, rradi);

%
% generate dictionary.
%

dict0 = [];

u_cordi0 = [];

dez0 = [];

xe0 = [];

ye0 = [];

[kr, kc] = size(x);


for i = 1:kr
    
    for j = 1:kc
        
        % e represents for 'element'
        [dict_e, e_size, xxe, yye] = gaussian_2Dspec(r_spec,s_spec,x(i,j),y(i,j));
        
        dict0 = [dict0 dict_e];
        
        u_cordi0 = [u_cordi0 [x(i,j) y(i,j)]'];
        
        dez0 = e_size;
        
        xe0 = xxe;
        
        ye0 = yye;
    end
    
end

dict = dict0; 

u_cordi = u_cordi0;

dez = dez0;

xe = xe0;

ye = ye0;

end
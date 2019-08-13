clear all;
clc;

% This is an dictionary generation example, the dictionary informations are
% saved to use later. Moreover, this file plots some examples of spectrum
% density representations.

[dict, dez, xe, ye] =  dic_gen(5,0.1,6,0.3);

columns = size(dict);

% plot spectrum density representations
for i = 1:2
    
    k = randperm(columns(1,2),1);
    
    figure
    pcolor(xe,ye,reshape(dict(:,k),dez))
    shading interp
    
end

% an structure saving variables
s.dict = dict;

s.dez = dez;

s.xe = xe;

s.ye = ye;


save('dict_infos.mat','-struct', 's');

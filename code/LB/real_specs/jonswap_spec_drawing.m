% This file simulates the directional wave spectrum modeled by the JONSWAP
% spectrum family
clear all
close all
clc

%% 
% simulation parameters
%
H = 2;% significant wave height(the average height of the 1/3 highest waves).(10)
        % this parameter don't change the shape of the frequency spectrum.
        % increase the energy of the spectrum.
       
T = 15;% the significant wave period.(14)
        % this parameter don't change the shape of the frequency spectrum.
        % move the spectrum towards the high frequencies.

c = 10;% Mitsuyasu-type spreading function's parameter.(10)

gam = 3.3;% the shape parameter of the JONSWAP frequency spectrum(3.3)
            % -- predicted parameter.

%x = 19.2639; % the fetch length (m)

%g = 9.81; % the gravitional constant (mps2)
%%
% retrieve the spectra's information
%
addpath(genpath('/homes/v18nguye/Documents/intern2019/data/WW3'))

str_date = '201001';
pnt_name = 'NODE008919';

% spectrum file
SPC_file = sprintf('ww3.%s_%s_spec.nc', pnt_name, str_date);
% image saved file name
sfn = sprintf('wind_sea_%s_%s.png',pnt_name,str_date);

% corresponding integrated wave parameters timeseries file
IWP_file = sprintf('ww3.%s_%s_para.nc', pnt_name, str_date);

% load files
SPC = Read_NetCDF(SPC_file);
IWP = Read_NetCDF(IWP_file);

% mask for found wave system(s)
mask_p0 = ~isnan(IWP.phs0); % partition 0 is for wind sea
mask_p1 = ~isnan(IWP.phs1); % partition 1 is for most energetic swell
mask_p2 = ~isnan(IWP.phs2); % partition 2 is for second most energetic swell
mask_p3 = ~isnan(IWP.phs3); % ....

% select date with detected wind sea, without swell
MatTime_WS = IWP.MatTime(mask_p0 & ~mask_p1 & ~mask_p2 & ~mask_p3);

% select date of first spectrum
MatTime = MatTime_WS(1);

% time indices in files
b1 = SPC.MatTime;
b2 = IWP.MatTime;
i1 = find(abs(SPC.MatTime-MatTime) < 1e-10);
i2 = find(abs(IWP.MatTime-MatTime) < 1e-10);

% get spectrum
d = SPC.direction([1:end,1]);
%freq  = SPC.frequency;
ang = SPC.direction([1:end,1]);
%theta = mod(-90-SPC.direction([1:end,1]),360) * pi/180;
Efth = SPC.efth([1:end,1],:,i1);
% the wave velocity
%Cp = SPC.cur(1,i1);
% the win velocity at 10 m elevation above the sea level (m/s)
%U = SPC.wnd(1,i1);
%U = 11;

freq = linspace(0,1.5,100);
theta = linspace(0,2*pi,36);
% creat a mesh grid for the frequency and the theta
[ffreq,ttheta] = meshgrid(freq,theta);

%% 
% simulation.
%
betaj = 0.06238*(1.094-0.01915*log(gam))/(0.23+0.0336*gam-0.185*((1.9+gam)^-1));

% nondimensional fetch
%xn = (g)*x/((U)^2);

% the peak of JONSWAP spectrum
%wp = 22*(g/U)*(nthroot(xn,-0.33));

% spectral peak
wp =2*pi*(1-0.132*((gam+0.2)^-0.559))/T;

% the Phillips constants
alpha = betaj*(H^2)*(wp^4);
%alpha = 0.076*(nthroot(xn,-0.22)); 

% jump function
w = ffreq;
sigma = (w/wp >1)*0.09 +(w/wp <= 1)*0.07;

% JONSWAP frequency spectrum function
phi = exp(-((w-wp).^2)./(2*(sigma.^2)*(wp^2)));
S_w = (alpha./(w.^5)).*exp(-1.25*((wp./w).^4)).*(gam.^phi);

%
% Mitsuyasu-type spreading function
%

% the wave age
%w_a = Cp/U;
%w_a = 0.52;

% spreading parameter
%smax = 11.5*(w_a^2.5);
%s = smax*((((w/wp > 1).*(w/wp)).^-2.5)+(((w/wp <= 1).*(w/wp)).^5));
%Go = (1/pi)*(2.^(2*s-1)).*(((gamma(s+1)).^2)./gamma(2*s+1));

% the mean wave direction
theta0 = pi/3;

fun = @(x) (cos((x- theta0)/2)).^(2*c);
Go = 1/(integral(fun, 0, 2*pi)); % normalizing parameter

%thetha0 = (SPC.curdir(1,i1)/180)*pi;
ttheta = mod(90-ttheta,360);
%s = 2;
G_theta = Go*((cos((ttheta-theta0)./2)).^(2*c));
% directional spectrum
S_w_theta = S_w.*G_theta;


%%
% plot the directional spectrum
%
figure('Name',sprintf('Directional Spectrum'))
[fx,fy] = pol2cart(ttheta,ffreq);
pcolor(fx,fy,S_w_theta)
%contour(fx,fy,S_w_theta)
shading flat
cb = colorbar;
set(get(cb,'ylabel'),'string','E(f,th) [m^2/Hz/rad]')

figure('Name',sprintf('JONSWAP freqency spectrum'))
plot(w(1,:),S_w(1,:))

figure('Name',sprintf('Mitsuyasu-type spreading function'))
plot(ttheta(:,1),G_theta(:,1))
function HIS = Read_datawell_his(datawell_his_file)
%
% Function Read_datawell_his.m reads the *.his file produced by W@ve21
% DATAWELL software and returns values in a structure.
%
% >> HIS = Read_datawell_his(datawell_his_file)
%
% returns the structure HIS  with fields :
%
%   MatDate : date in matlab format [days]
%   Tp      : the peak period (the reciprocal of the peak frequency) [s]
%   Dirp    : the wave direction at the peak frequency [°]
%   Sprp    : the directional spread at the peak frequency [°]
%   Tz      : the zero-upcross period [s]
%   Hs      :  the significant wave height [cm]
%   TI      : the integral period, or Tm(-2,0) [s]
%   T1      : the mean period, or Tm(0,1) [s]
%   Tc      : the crest period, or Tm(2,4) [s]
%   Tdw2    : wave period Tm(-1,1) [s]
%   Tdw1    : peak period estimator [s]
%   Tpc     : calculated peak period [s]
%   nu      : Longuet-Higgins bandwidth parameter [-]
%   eps     : bandwidth parameter [-]
%   QP      : Goda’s peakedness parameter [-]
%   Ss      : significant steepness [-]
%   Tref    : reference temperature []
%   Tsea    : Sea surface temperature []
%   Bat     : battery status (0 = empty to 7 = full)
%

% check input file
if ~exist(datawell_his_file,'file')
    error('Read_datawell_his: File %s not found', datawell_his_file)
end

% check input file name
[~,~,ext] = fileparts(datawell_his_file);
if ~strcmp(ext,'.his')
    warning(['Extension of input file is not the expected one (.his) \n',...
        'Program may result in error or bad data reading'])     
end

% open file
fid = fopen(datawell_his_file,'r');

% format
Nvals_per_line = 24;
format = '%d-%d-%dT%d:%d:%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d';


% read file
[data, Nvals] = fscanf(fid,format,Inf);

% close file
fclose(fid);

% reshape
Nlines = Nvals / Nvals_per_line;
data = reshape(data,Nvals_per_line,Nlines);

% get values and store them in output structure
HIS.MatDate = datenum(data(1:6,:)');
HIS.Tp      = data( 7,:)';
HIS.Dirp    = data( 8,:)';
HIS.Sprp    = data( 9,:)';
HIS.Tz      = data(10,:)';
HIS.Hs      = data(11,:)';
HIS.TI      = data(12,:)';
HIS.T1      = data(13,:)';
HIS.Tc      = data(14,:)';
HIS.Tdw1    = data(15,:)';
HIS.Tdw2    = data(16,:)';
HIS.Tpc     = data(17,:)';
HIS.nu      = data(18,:)';
HIS.eps     = data(19,:)';
HIS.QP      = data(20,:)';
HIS.Ss      = data(21,:)';
HIS.Tref    = data(22,:)';
HIS.Tsea    = data(23,:)';
HIS.Bat     = data(24,:)';

return

fprintf('MatTime : date in matlab format [days]\n');
fprintf('Tp      : the peak period (the reciprocal of the peak frequency) [s]\n');
fprintf('Dirp    : the wave direction at the peak frequency [°]\n');
fprintf('Sprp    : the directional spread at the peak frequency [°]\n');
fprintf('Tz      : the zero-upcross period [s]\n');
fprintf('Hs      :  the significant wave height [cm]\n');
fprintf('TI      : the integral period, or Tm(-2,0) [s]\n');
fprintf('T1      : the mean period, or Tm(0,1) [s]\n');
fprintf('Tc      : the crest period, or Tm(2,4) [s]\n');
fprintf('Tdw2    : wave period Tm(-1,1) [s]\n');
fprintf('Tdw1    : peak period estimator [s]\n');
fprintf('Tpc     : calculated peak period [s]\n');
fprintf('nu      : Longuet-Higgins bandwidth parameter [-]\n');
fprintf('eps     : bandwidth parameter [-]\n');
fprintf('QP      : Goda’s peakedness parameter [-]\n');
fprintf('Ss      : significant steepness [-]\n');
fprintf('Tref    : reference temperature []\n');
fprintf('Tsea    : Sea surface temperature []\n');
fprintf('Bat     : battery status (0 = empty to 7 = full)\n');

end


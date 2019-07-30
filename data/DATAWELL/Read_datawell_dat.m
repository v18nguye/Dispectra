function DAT = Read_datawell_dat(datawell_dat_file)
%
% Function Read_datawell_dat.m reads the *.dat file produced by W@ve21
% DATAWELL software and returns read spectrum.
%
% >> DAT = Read_datawell_dat(datawell_dat_file)
%
%  returns the structure DAT  with fields :
%
%    freq : frequencies [Hz]
%    dir  : direction [deg] 
%    Sfth : spectral density [m^2/Hz/rad]
% 
%    MatTime : Date (matlab format) extracted from file name [only if 
%              W@ve21 format is respected] 
%

% check input file
if ~exist(datawell_dat_file,'file')
    error('Read_datawell_dat: File %s not found', datawell_dat_file)
end

% check input file name
[unused,name,ext] = fileparts(datawell_dat_file);
if ~strcmp(ext,'.dat')
    warning(['Extension of input file is not the expected one (.dat) \n',...
        'Program may result in error or bad data reading'])     
end

% try to extract date from file name (ex: 2015-10-22T09h23.dat)
try 
    [V,n] = sscanf(name,'%d-%d-%dT%dh%d');
    if n ~= 5; error('Bad filename'); end
    MatDate0 = datenum(V(1),V(2),V(3),V(4),V(5),0);
    %fprintf('File date : %s \n', datestr(MatDate0))
catch 
    MatDate0 = NaN;
    warning(['Program can extract date from input filename :\n',...
        'File name does not not respect W@ve21 format']) 
end

% open file
fid = fopen(datawell_dat_file);

% read header
line = fgets(fid);
C = textscan(line,'%s %d %d');

% check file format
if ~strcmp(C{1}{1},'IRGRID')
    error('Read_datawell_dat: Bad file format')
end

% spectrum size
Nfreq = double(C{2});
Ndir  = double(C{3});

% read header
fgets(fid);

% read frequencies
DAT.freq = fscanf(fid,'%f',Nfreq)';

% read directions
DAT.dir = fscanf(fid,'%f',Ndir);

% read spectrum
E10 = NaN(Ndir,Nfreq);
for ifreq = 1:Nfreq
    E10(:,ifreq) = fscanf(fid,'%f',Ndir);
end
DAT.Sfth = 10.^E10;

% mask bad values
mask = E10 == -4;
DAT.Sfth(mask) = 0;

% add date
DAT.MatDate = MatDate0;

% close file
fclose(fid);

return

%% DEBUG

freq = DAT.freq;
dir  = DAT.dir;
Sfth = DAT.Sfth;

figure(999)
pcolor(freq,dir,Sfth)
xlabel('freq [Hz]')
ylabel('dir [deg]')
cb = colorbar;
set(get(cb,'xlabel'),'string',{'S','[m^2/Hz/rad]'})

end
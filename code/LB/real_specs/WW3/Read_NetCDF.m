function DATA = Read_NetCDF(NCfile, varargin)
%
% Fabien Leckler, fabien.leckler@shom.fr
%
% This function extracts data from NetCDF file and returned them in a
% structure. 
%
% The function processes extraction applying 'scale_factor', 'add_offset' 
% and '_FillValue' if those attributes are present in NetCDF file. 
%
% Moreover, the function try to return 'time' variable in Mat-Format 
% (i.e. datenum format) using the attribute 'units' of variable 'time' :
% if time units format is recognized, the function return an additionnal 
% field MatTime.
%
% Please report bugs or ask for improvement to fabien.leckler@shom.fr 
%
% USAGE : data = Read_NetCDF(ncid, [OPTION])
% 
% INPUT:
% ------
%   ncid or filename
%          NetCDF file identifier returned by ncid = netcdf.open(filename). 
%          However, fuction also accept in first argument the 'string' name
%          of NetCDF file.
%   
% OUTPUT:
% -------
%   DATA   Structure
%
% OPTION:
% -------
%   '-dim', STRDIM1, [F1:D1:E1], STRDIM2, [F2:D2:E2], ... 
%
%             This option is used to select part of variable to extract
%             - STRDIM  is the dimension name.
%             - [F:D:E] is the indice(s) to extract. WARNING: Vector 
%                       must be monotonic and ascending;
%
%   '-var', STRVAR1, STRVAR2, ...
%
%            This option is used to select variable to extract, default
%            extract all variables.
%          
%
% EXEMPLES: 
% 
% Here, File.nc contains 2 dimensions X and Y, and 4 variables V1(X),
% V2(Y), V3(X,Y) and V4(X,Y).
%
% >> data = Read_NetCDF('File.nc')
% data =
%    v1: V1(1:end)
%    v2: V2(1:end)
%    v3: V3(1:end,1:end)
%    v4: V4(1:end,1:end)
%
% >> data = Read_NetCDF('File.nc','-var', 'V1','V2') 
% data =
%    v1: V1(1:end)
%    v2: V2(1:end)
%
% >> ncid = netcdf.open('File.nc');
% >> data = Read_NetCDF(ncid,'-var', 'V1','V2') 
% data =
%    v1: V1(1:end)
%    v2: V2(1:end)
%
% >> data = Read_NetCDF('File.nc','-dim','X', 2, 'Y', 3) 
%  data =
%    v1: V1(2)
%    v2: V2(3)
%    v3: V3(2,3)
%    v4: V4(2,3)
%
% >> data = Read_NetCDF('File.nc','-dim','X',[2:10],'-var','DATA') 
%  data =
%    v1: V1(2:10)
%    v2: V2(1:end)
%    v3: V3(2:10,1:end)
%    v4: V4(2:10,1:end)
%
% >> data = netcdf.getVar('File.nc','-var', 'V1', 'V3','-dim','Y', [2:3:31]) 
%  data =
%    v1: V1(1:end)
%    v3: V3(1:end,2:3:31)
%
%

% get file id or file name
if isnumeric(NCfile)
    ncid  = NCfile;
    CLOSE = false; 
elseif ischar(NCfile) 
    if exist(NCfile,'file')
        ncid = netcdf.open ( NCfile,'NOWRITE');
        CLOSE = true; 
    else
        error(['Read_NetCDF: Bad Input argument (1):' ...
            'File %s not found'], NCfile)
    end
else
    error(['Read_NetCDF: Bad Input argument (1):' ...
        'File ID or File Name expected'])
end

% Read file and get number of dimensions and of variables
%%% netcdf.inq returns information [ndims,nvars,ngatts,unlimdimid] = number
%%% of dimensions, variables.
[NumDim, NumVar, ~, ~] = netcdf.inq(ncid);
for iDim=1:NumDim
    %%% netcdf.inqDim returns name and length of dimensions specified by
    %%% iDim -1.
    [DimName{iDim}, DimLen(iDim)] = netcdf.inqDim(ncid,iDim-1);
end
for iVar=1:NumVar
    %%% return information of a variable specified by ivar - 1 = name, ~,
    %%% dimension IDs, and the number of attributes associated with variable. 
    [VarName{iVar},~,Dimids{iVar},Natts{iVar}] = netcdf.inqVar(ncid,iVar-1);
end

NumDim; % NumDim = 5
DimLen; % DimLen = 744     1    16    30    24

VarName; % {'time'}    {'station'}    {'string16'}    {'station_name'}    {'longitude'}    {'latitude'}    {'frequency'}    {'frequency1'}    {'frequency2'}
% {'direction'}    {'efth'}    {'dpt'}    {'wnd'}    {'wnddir'}    {'cur'}    {'curdir'}
Dimids; % each varName have a diffrent dimension.
%     {[0]}    {[1]}    {[2]}    {1×2 double}    {1×2 double}    {1×2 double}    {[3]}    {[3]}    {[3]}    {[4]}    {1×4 double}    {1×2 double}
% {1×2 double}    {1×2 double}    {1×2 double}    {1×2 double}
Natts; % {[6]}    {[6]}    {[5]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}


% Default : all values are extract (converted for netcdf.getVar)
INDEX  = zeros([1,NumDim]);
LENGTH = DimLen;
STRIDE = ones([1,NumDim]);

% Initialize variables to extract
VarMask = true([1,NumVar]);

% Number of varargin arguments
%%% varagin is left for optinal arguments
Narg = numel(varargin);

% Read Option 
if Narg > 0

    % first varargin argument must be an option 
    if ~ischar(varargin{1}) || ~strcmp(varargin{1}(1),'-')
        error(['Read_NetCDF: Bad Input argument (2):' ...
            'Option argument expected'])
    end
    
    % variable defined by user ?
    if any(strcmp(varargin{1},'-var'))
        VarMask(:) = false;
    end

    % loop over argument
    iarg = 1;
    while iarg <= Narg

        if strcmp(varargin{iarg},'-var')

            % read variable name
            while iarg+1 <= Narg && ~strcmp(varargin{iarg+1}(1),'-')
                
                % next argument 
                iarg = iarg+1;

                % check if string
                if ~ischar(varargin{iarg}) 
                    error(['Read_NetCDF: Bad Input argument (%d):' ...
                        'String argument expected'], iarg+1)
                end
                
                % get variable
                if any(strcmp(varargin{iarg},VarName))
                    iVar = netcdf.inqVarID(ncid,varargin{iarg})+1;
                    VarMask(iVar) = true;
                else
                    error(['Read_NetCDF: Bad argument (%d):' ...
                        'Variable not in file (%s)'], iarg+1, varargin{iarg});
                end    

            end
            
        elseif strcmp(varargin{iarg},'-dim')

            % read dimention custom
            while iarg+1 <= Narg && ~strcmp(varargin{iarg+1}(1),'-')
                
                % next argument 
                iarg = iarg+1;

                % check if string
                if ~ischar(varargin{iarg}) 
                    error(['Read_NetCDF: Bad Input argument (%d):' ...
                        'String argument expected'], iarg)
                end
                
                % get dim
                if any(strcmp(varargin{iarg},DimName))
                    iDim = netcdf.inqDimID(ncid,varargin{iarg})+1;
                else
                    error(['Read_NetCDF: Bad argument (%d):' ...
                        'Dimension not in file (%s)'], iarg, varargin{iarg});
                end    
                
                % next argument
                if iarg+1 <= Narg
                    iarg = iarg+1;
                else
                    error(['Read_NetCDF: Missing argument :' ...
                        'Vector containing index to extract for ' ...
                        'dimention %s is expected'], varargin{iarg});
                end    
                
                % check if numeric
                if ~isnumeric(varargin{iarg}) || ~isvector(varargin{iarg}) 
                    error(['Read_NetCDF: Bad Input argument (%d):' ...
                        'Numeric vector argument expected'], iarg)
                end
                
                % get index(s) to extract
                Vec = varargin{iarg};
                if numel(Vec) == 1
                
                    % convert index to  netcdf.getVar argument 
                    INDEX (iDim) = Vec-1;
                    LENGTH(iDim) = 1;
                    STRIDE(iDim) = 1;
                    
                else

                    % check if monotonic and ascending
                    D = mean(diff(Vec));
                    if any(diff(Vec)~=D) || D < 0
                        error(['Read_NetCDF: Bad Input argument (%d):' ...
                            'Indice vector to extract must be ' ...
                            'monotonic and ascending'], iarg)
                    elseif D ~= round(D)
                        error(['Read_NetCDF: Bad Input argument (%d):' ...
                            'Indice to extract must be integer'], iarg)
                    end
                    
                    % convert index to  netcdf.getVar argument 
                    INDEX (iDim) = Vec(1)-1;
                    LENGTH(iDim) = numel(Vec);
                    STRIDE(iDim) = D;
                    
                end
            end
            
        else
            
            error(['Read_NetCDF: Bad Input argument (%d):' ...
                'Not Reconized'], iarg)
            
        end
        
        % next argument
        iarg = iarg+1;

    end
end

% Any variable to extract ?
if ~any(VarMask)
    error(['Read_NetCDF: Bad Input argument:' ...
        'No variable to output'])
end

% Extract variable(s)
%%% for general case
for iVar = 1:NumVar
    
    % variable to extract ?
    %%% if there are no variable at VarMask(ivar), pass to the next loop
    %%% control.
    if ~VarMask(iVar)
        continue
    end
    
    % specific case for time (out time in matlab-date)
    %%% strcmpi: compare 2 strings
    if strcmpi(VarName{iVar}, 'time')
        DATE = true;
        sf2  = NaN;
        ao2  = NaN;
    else
        DATE = false;
    end
    VarName; % {'time'}    {'station'}    {'string16'}    {'station_name'}    {'longitude'}    {'latitude'}    {'frequency'}    {'frequency1'}    {'frequency2'}
% {'direction'}    {'efth'}    {'dpt'}    {'wnd'}    {'wnddir'}    {'cur'}    {'curdir'}
    Dimids; % each varName have a diffrent dimension.
%     {[0]}    {[1]}    {[2]}    {1×2 double}    {1×2 double}    {1×2 double}    {[3]}    {[3]}    {[3]}    {[4]}    {1×4 double}    {1×2 double}
% {1×2 double}    {1×2 double}    {1×2 double}    {1×2 double}

    Natts; % {[6]}    {[6]}    {[5]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}    {[9]}

    
    % dimension of var
    %%% Dimids
    iDims = Dimids{iVar}+1;
    start  = INDEX (iDims);
    count  = LENGTH(iDims);
    stride = STRIDE(iDims);
    
    % get variable 
    data = double(netcdf.getVar(ncid,iVar-1,start,count,stride));
    
    % default
    sf = 1;
    ao = 0;
    fv = NaN;
    
    %loop over attribute
    for iAtt = 1:Natts{iVar}
        % return the attname
        attname = netcdf.inqAttName(ncid,iVar-1, iAtt-1);
        
        if strcmp(attname, 'scale_factor')
            sf    = sf * double(netcdf.getAtt(ncid, iVar-1, 'scale_factor'));
        
        elseif strcmp(attname, 'add_offset')
            ao    = ao + double(netcdf.getAtt(ncid, iVar-1, 'add_offset'));
        
        elseif strcmp(attname, '_FillValue')
            fv    = double(netcdf.getAtt(ncid, iVar-1, '_FillValue'));
        
        elseif DATE && strcmp(attname, 'units')
            
            % get units
            units = netcdf.getAtt(ncid, iVar-1, 'units');
            if ~isempty(strfind(lower(units),'second'))
                sf2 = 1 / (24*3600);
            elseif strfind(lower(units),'minutes')
                sf2 = 1 / (24*60);
            elseif strfind(lower(units),'hours')
                sf2 = 1 / (24*1);
            elseif strfind(lower(units),'days')
                sf2 = 1 ;
            end
            
            % get reference date
            if strfind(units,'since')
                refdate =  units(strfind(units,'since')+5:end);
                try
                    ao2 = datenum(refdate);
                catch me
                    try
                        % for WW3
                        refdate2 = refdate;
                        refdate2(refdate2 == 'Z') = ' ';
                        refdate2(refdate2 == 'T') = ' ';
                        ao2 = datenum(refdate2);
                    end
                end
            end
            
        end
    end
    
    % apply NetCDF scale facor and add offset
    mask  = data == fv;
    data  = ao + double(data) * sf;
    data(mask) = NaN;
    
    eval(['DATA.' VarName{iVar} '= squeeze(data);']);

    % apply scale facor and add offset from unit for date
    if DATE && ~isnan(ao2) && ~isnan(sf2)
        DATA.MatTime  = ao2 + double(data) * sf2;
    end

end

if CLOSE
    netcdf.close(ncid)
end

end

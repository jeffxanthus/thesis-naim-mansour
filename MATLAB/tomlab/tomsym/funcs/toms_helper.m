function y = toms_helper(fname, args, varargin)
% toms_helper - Helper function for toms
%
% This function is not meant to be called directly.

% y = toms_helper(fname, args, varargin) creates a cell array y so that y{i}
% corresponds to the name in varargin{i}, if that is not a flag or size
% specification. The function fname is used to create the objects.
%
% This function is called by toms, tomStates and tomControls

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-09-08 by rutquist for TOMLAB release 7.7

y = cell(size(varargin));

sz1 = 1;
sz2 = 1;

flags = {};

inds = false;

isint = false;
issym = false;

for i=1:length(varargin)
    % Parse for strings like "2x3", "int", "symmetric"
    vi = varargin{i};
    if any(strcmp(vi,{'-int','-integer','int','integer'}))
        flags = {flags{:},'int'};
        isint = true;
        continue
    end
    if any(strcmp(vi,{'-symmetric','symmetric'}))
        flags = {flags{:},'symmetric'};
        issym = true;
        continue
    end
    xi = find(vi=='x');
    if length(xi)==1
        szt1 = str2double(vi(1:xi-1));
        if vi(end)=='!'
            ig = true;
            szt2 = str2double(vi(xi+1:end-1));
        else
            ig = false;
            szt2 = str2double(vi(xi+1:end));
        end
        if ~any(isnan([szt1 szt2]))
            sz1 = szt1;
            sz2 = szt2;
            inds = ig;
            continue;
        end
    end

    if ~isvarname(varargin{i})
        error(['Illegal variable name ' varargin{i}]);
    end
    
    if inds && sz1*sz2 > 1
        % Create a matrix from scalar symbols
        if sz1>9 || sz2>9
            div = '_';
        else
            div = '';
        end
        if issym
            if isint
                flg = {'int'};
            else
                flg = {};
            end
            if sz1~=sz2
                error('Symmetric matrix must be square');
            end
            isy = cell(1,sz1*(sz1+1)/2);
            [iix,iiy] = meshgrid(1:sz1,1:sz2);
            eli = find(tril(ones(sz1)));
            for ii=1:length(isy);
                isy{ii} = feval(fname,args{:},[varargin{i} div num2str(iix(eli(ii))) div num2str(iiy(eli(ii)))],1,1,flg{:});
            end
            isv = setSymmetric(vertcat(isy{:}));
        elseif sz1==1 || sz2==1
            isy = cell(1,sz1*sz2);
            for ii=1:length(isy);
                isy{ii} = feval(fname,args{:},[varargin{i} num2str(ii)],1,1,flags{:});
            end
            if sz1>0
                isv = vertcat(isy{:});
            else
                isv = horzcat(isy{:});
            end
        else
            isy = cell(1,sz1);
            for ii=1:length(isy);
                isx = cell(1,sz2);
                for j=1:length(isx)
                    isx{j} = feval(fname,args{:},[varargin{i} div num2str(j) div num2str(ii)],1,1,flags{:});
                end
                isy{ii} = vertcat(isx{:});
            end
            isv = horzcat(isy{:});
        end
        y{i} = isv;
    else
        y{i} = feval(fname,args{:},varargin{i},sz1,sz2,flags{:});
    end
end

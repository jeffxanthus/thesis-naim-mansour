function y = wrap(info,varargin)
% tomSym/wrap - Overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-05-17 by rutquist for TOMLAB release 7.7

% Fix common capitalization mistake
if isfield(info,'JFuns') && ~isfield(info,'Jfuns')
    info.Jfuns = info.JFuns;
    info = rmfield(info,'JFuns');
end
if ~isfield(info,'Jpatterns') && isfield(info,'JPatterns')
    info.Jpatterns = info.JPatterns;
    info = rmfield(info,'JPatterns');
end

if isfield(info,'Jpatterns')
    info.Jcolidx = cell(size(info.Jpatterns));
    for i=1:length(info.Jpatterns)
        info.Jcolidx{i}   = findpatt(info.Jpatterns{i});
    end
end

% Expand single string into a list.
if isfield(info,'Jfuns') && ischar(info.Jfuns)
    str = info.Jfuns;
    info.Jfuns = cell(1,length(varargin));
    for i=1:length(info.Jfuns)
        info.Jfuns{i} = str;
    end
end

% Fill up with empty if no Jpatterns given
if ~isfield(info,'Jpatterns')
    info.Jpatterns = cell(1,length(varargin));
    for i=1:length(varargin)
        info.Jpatterns{i} = [];
    end
end

% Fill up with empty if no Jcolidx given
if ~isfield(info,'Jcolidx')
    info.Jcolidx = cell(1,length(varargin));
    for i=1:length(varargin)
        info.Jcolidx{i} = [];
    end
end
    
y = tomSym(mfilename,info.sz1,info.sz2,info,varargin{:});

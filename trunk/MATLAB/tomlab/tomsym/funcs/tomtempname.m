function name = tomtempname(flag)
% tomtempname - Get a temporary file name.
%
% tomtempname returns a unique name for use as a tomSym temporary file.
%
% tomtempname depends on Matlab's native TEMPNAME, but modifies the
% filename to contain a marker, so that the files can be identified for
% cleanup.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

marker = 'tomSymTmp';

if nargin>0
    if ischar(flag) && strcmpi(flag,'marker');
        name = marker;
        return
    else
        error('Unrecognized flag.');
    end
end 

ok = false;
for i=1:3
    name = [tempname num2str(floor(10*rand)) marker];
    if ~exist(name,'file');
        ok = true;
        break
    end
end

if ~ok
    error('Failed to generate a unique file name.');
end

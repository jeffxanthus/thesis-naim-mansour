function X = double(M)
% tomArray/double - convert a tomArray to double
%
%  X = double(M) converts tomArray M into an ordinary Matlab array.
%
% See also: tomArray/unArray

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~isnumeric(M.X)
    error('Cannot convert tomArray into double. Data is not numeric.');
end

X = unArray(M);

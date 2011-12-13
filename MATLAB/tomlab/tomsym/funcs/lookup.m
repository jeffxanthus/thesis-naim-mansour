function y = lookup(v, i)
% Lookup - Extract elements from a matrix
%
% y = lookup(v, i) is equivalent to y = v(i).
%
% lookup(v,':') is an alias for vec(v)
%
% This function is defined so that it can be overloaded for tomSym,
% since Matlab syntax rules make some constructions impossible otherwise.
%
% See also submatrix

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-14 by rutquist for TOMLAB release 7.7

if ischar(i) && strcmp(i,':')
    y = vec(v);
else
    y = v(i);
end

function X = unArray(M)
% unArray - extract data from a tomArray
%
%  X = unArray(M) converts tomArray M into an ordinary Matlab array or a
%  tomSym symbolic expression.
%
% Numeric tomArrays of any size can be converted into Matlab arrays,
% keeping their size.
%
% When the tomArray contains tomSym data, it is only possible to use
% unArray if the number of dimensions is one or two. Arrays of higher
% dimension must first be reudced using, for example, vec, reshape or
% squeeze.
%
% See also: tomArray/vec

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if isa(M.X,'tomSym') && length(M.sz)>2
    error('Too many dimensions for conversion to tomsym.');
end

if isempty(M.sz)
    M.sz = [1 1];
elseif length(M.sz)==1
    M.sz = [M.sz 1];
end

X = reshape(M.X, M.sz);

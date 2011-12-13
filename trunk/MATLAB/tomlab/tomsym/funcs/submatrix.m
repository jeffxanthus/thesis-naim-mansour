function y = submatrix(M, i, j)
% Lookup - Extract rows/columns from a matrix
%
% y = submatrix(M, i, j) is equivalent to y = M(i,j)
%
% i and j must be vectors of indexes or the character ':'
% (Note that a colon must be passed using quote marks, 
% for example lookup(M, 1, ':'))
%
% This function is defined so that it can be overloaded for tomSym.
%
% See also lookup

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-14 by rutquist for TOMLAB release 7.7

if ischar(i) && strcmp(i,':')
    i = 1:size(M,1);
end
if ischar(j) && strcmp(j,':')
    j = 1:size(M,2);
end
ri = round(i);
rj = round(j);
if max(abs(ri-i)) > 1e-9 || max(abs(rj-j))
    error('Subscript indices must either be real positive integers or logicals.');
end

y = M(ri,rj);

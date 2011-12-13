function z = tomzeros(varargin)
% tomzeros - Short for tomSym(zeros(...))
%
% z = tomzeros(n) gives an n-by-n matrix of zeros.
%
% z = tomzeros(m,n) or tomzeros([m,n]) is an m-by-n matrix of zeros.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

z = tomSym(zeros(varargin{:}));

function y = linspace(x1, x2, n)
% tomSym/linspace - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<3
    n = 100;
end

if isa(n,'tomSym')
    error('N must be a numeric value.');
end

y = x1+(x2-x1)*linspace(0,1,n);



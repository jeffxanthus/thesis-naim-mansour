function y=spower(a,b)
% spower - Element-wise power to a scalar exponent.
%
% Y=SPOWER(A,B) is the same as Y=A.^B if A is a matrix and B is a scalar.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(b) ~= 1
    error('Spower expects a scalar exponent.');
end

y = a.^b;

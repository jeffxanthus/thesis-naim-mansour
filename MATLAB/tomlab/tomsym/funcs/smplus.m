function y=smplus(a,b)
% smtimes - Scalar + matrix addition
%
% This operation maps to plus, and uses that operator, but in tomSym it
% is treated separate, because it is not the same thing as matrix+matrix
% addition.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = a+b;

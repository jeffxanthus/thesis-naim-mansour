function y = nOperands(o)
% tomSym/nOperands - Get number of operands from a tomSym
%
% y = nOperands(o) gets the number of operands from tomSym object o.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = length(o.s(end).a);

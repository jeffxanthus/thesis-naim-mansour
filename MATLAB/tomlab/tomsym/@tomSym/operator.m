function y = operator(o)
% tomSym/operator - Get the operator from a tomSym
%
% y = operand(o) gets the operator text string from tomSym object o.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

y = o.s(end).op;

function y=smtimes(a,b)
% smtimes - Scalar x matrix multiplication
%
% This operation maps to mtimes, and uses that operator, but in tomSym it
% has to be treated separate, because it is not the same thing as matrix
% multiplication (It is the same thing as kron, but is kept separate from
% that too, in order to be displayed as in matlab). The scalar is always
% moved to the left.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = a*b;

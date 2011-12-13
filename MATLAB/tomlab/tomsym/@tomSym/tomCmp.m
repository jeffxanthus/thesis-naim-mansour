function y = tomCmp(o, s)
% tomCmp - Test if the operator of a tomSym is of a certain type
%
% Example: tomCmp(x+2, 'plus') returns true.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-28 by rutquist for TOMLAB release 7.7

% Overloaded version of the function is called if o is a tomSym.

y = strcmp(o.s(end).op,s);

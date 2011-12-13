function y = tomCmp(o, s)
% tomCmp - Test if the operator of a tomSym is of a certain type
%
% Example: tomCmp(x+2, 'plus') returns true.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-05 by rutquist for TOMLAB release 7.7

% If o is really a tomSym, then the overloaded version of this function
% will be called.

y = false;

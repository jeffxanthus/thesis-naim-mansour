function y = tomCmpO(n, o, s)
% tomCmp - Test if the operator of a tomSym operand is of a certain type
%
% tomCmp(n,o,s) does the same thing as tomCmp(operand(n,o),s) but faster.
%
% Example: tomCmpO(2,3+cos(x), 'cos') returns true.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% $Id $

% Overloaded version of the function is called if o is a tomSym.

y = o.s(end).a(n) > 0 && strcmp(o.s(o.s(end).a(n)).op,s);

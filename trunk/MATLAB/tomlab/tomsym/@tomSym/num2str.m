function c = num2str(p)
% tomSym/num2str - Convert tomSym to a string.
%
% A tomSym is not, strictly speaking, a number. But overloading this
% function does what the user expects.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

c = char(p);

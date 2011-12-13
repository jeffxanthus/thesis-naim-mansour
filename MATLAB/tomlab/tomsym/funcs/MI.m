function y = MI(a) %#ok
% MI - Matrix Inequality (Only for tomSym objects)
%
% MI cannot be used with non-tomSym matrices, because the <= and >=
% operators already have a different meaning for numeric matrices.
%
% See: tomSym/MI

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

error('MI cannot be used on a constant.');

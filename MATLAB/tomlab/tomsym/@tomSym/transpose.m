function p = transpose(a)
% tomSym/transpose - Overloaded operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

% TomSym makes no distinction between ctranspose and transpose, because
% it does not handle complex numbers.

p = ctranspose(a);

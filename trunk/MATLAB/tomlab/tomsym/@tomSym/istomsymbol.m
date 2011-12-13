function y = istomsymbol(a)
% tomSym/istomsymbol - Returns "true" if a is a tomSym symbol.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = strcmp(operator(a),'tom');

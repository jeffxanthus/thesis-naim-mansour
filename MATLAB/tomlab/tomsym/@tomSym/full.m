function a = full(a)
% tomSym/full - Overloaded function
%
% If the tomSym is a consant, then the result is a full matrix.
%
% Otherwise "full" has no effect on a tomSym. The result may still be a
% sparse matrix when values are substituted for symbols.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if tomCmp(a,'constant')
    a = full(double(a));
end

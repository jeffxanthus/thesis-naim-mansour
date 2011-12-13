function y = eye(a)
% tomSym/eye - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(a)~=1
    error('Symbolic identity matrix only takes one input argument');
end

% Simplifications
if isone(a)
    y = 1;
else
    y = tomSym(mfilename, a, a);
end

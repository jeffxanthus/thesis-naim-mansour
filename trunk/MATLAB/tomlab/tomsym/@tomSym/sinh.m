function y = sinh(a)
% tomSym/sinh - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

% Simplifications
if strcmp(operator(a),'asinh')
    y = operand(1,a);
else
    y = quickop(mfilename,a);
end

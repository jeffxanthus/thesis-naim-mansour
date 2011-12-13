function y = acos(a)
% tomSym/acos - Overloaded function
%
% The function assumes that abs(a) < 1, so acos(cos(x)) == x

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Simplifications - We assume that a is constrained to the proper
% intervall.
if strcmp(operator(a),'cos')
    % Assuming that abs(a)<1.
    y = operand(1,a);
else
    y = quickop(mfilename,a);
end

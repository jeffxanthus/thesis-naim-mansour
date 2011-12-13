% Simple test function for fminunc, example page 4-61-4-62 in OPTB 2.0 manual
%
% function [f,g]=fminunc_fg2(x)
%
% f = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2
%
% g = [6*x(1) + 2*x(2); 2*x(1) + 2*x(2)];

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 29, 1999.   Last modified July 29, 1999.

function [f, g] = fminunc_fg2(x)

f = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2;
if nargout > 1
   g = [6*x(1) + 2*x(2); 2*x(1) + 2*x(2)];
end
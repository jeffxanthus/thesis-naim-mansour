function y = DiracDeltas(x,breaks,coefs,orders)
% tomSym/DiracDeltas - Placeholder derivative of a function with multiple discontinuities.
%
% y = DiracDeltas(x,breaks) represents the sum over i of DiracDelta(x-break(i)) 
%
% y = DiracDeltas(x) without specifying the position of the breaks
% represents a function that is zero everywhere but is not the derivative
% of a constant function. (TomSym uses this internally to diagnose
% optimization problems that would benefit from restructuring.)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    y = tomSym(mfilename,size(x,1),size(x,2),x);
else
    y = tomSym(mfilename,size(x,1),size(x,2),x,breaks,coefs,orders);
end

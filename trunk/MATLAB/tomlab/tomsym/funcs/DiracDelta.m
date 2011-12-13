function y = DiracDelta(x,order)
% DiracDelta - Placeholder derivative of a discontinuous function.
%
% The derivative of a discontinuous function, such as sign(x) is not
% well-defined in the point where the discontinuity occurs. In order to
% reflect this, tomSym uses the convention that
%    derivative(sign(x),x) = 2*DiracDelta(x)
% where DiracDelta(x) is a function that is zero everywhere, except at
% x=0.
%
% By doing so, tomSym can more easily diagnose and categorize problems that
% contain discontinuities.
%
% The derivative of DiracDelta(x) is DiracDelta(x,2), and in general the
% derivative of DiracDelta(x,n) is DiracDelta(x,n+1). Antiderivatives are
% also defined. The Heaviside function can be written as DiracDelta(x,0),
% and a ramp as DiracDelta(x,-1).
%
% See also: DiracDeltas

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if nargin<2
    order = 1;
end

y = spzeros(size(x));

if order<=0
    y(x>0) = x(x>0).^(-order) / factorial(-order);
else
    %y(x==0) = NaN; %This would be mathematically correct, but just
                    %troublesome when using numeric solvers.
end

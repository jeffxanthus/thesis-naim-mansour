function y = DiracDeltas(x,breaks,coefs,orders)
% DiracDeltas - Placeholder derivative of a function with multiple discontinuities.
%
% y = DiracDeltas(x,breaks,orders) represents the sum over i of 
% reshape(coefs(i,:),size(x)).*DiracDelta(x-breaks(i),orders(i))
%
% y = DiracDeltas(x) simply returns zero. This notation is sometimes used 
% by tomSym for diagnostic purpouses, to indicate the derivative of a 
% piecewise constant function.
%
% See also: DiracDelta

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = spzeros(size(x));
if nargin>2
    for i=1:numel(breaks)
        y = y + reshape(coefs(i,:),size(x)).*DiracDelta(x-breaks(i),orders(i));
    end
end

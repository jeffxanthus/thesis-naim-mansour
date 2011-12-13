function D=setdiag(v)
% tomSym/setdiag - Setdiag for tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(v)==1
    % Simplify if the argument is zero or scalar.
    D = v;
elseif tomCmp(v,'vec')
    D = setdiag(operand(1,v));
else
    D = tomSym(mfilename,numel(v),numel(v),v);
end

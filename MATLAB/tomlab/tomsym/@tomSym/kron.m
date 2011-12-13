function y=kron(a,b)
% tomSym/kron - Overload the Kronecker product

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-03 by rutquist for TOMLAB release 7.7

sz1 = size(a,1)*size(b,1);
sz2 = size(a,2)*size(b,2);

if iszero(a) || iszero(b)
    y = zeros(sz1,sz2);
elseif numel(a)==1 || numel(b)==1
    y = a*b;  % Will map to smtimes
elseif isone(a)
    y = ekron(b,size(a,1),1);
elseif isone(b)
    y = ekron(a,size(b,1),2);
else
    y = tomSym(mfilename,sz1,sz2,a,b);
end

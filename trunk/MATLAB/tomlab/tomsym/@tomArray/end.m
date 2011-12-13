function e = end(o,k,n)
% tomArray/end
% Overloaded - Get the last index of a tomArray object.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if n==1
    e = prod(o.sz);
else
    e = o.sz(k);
end

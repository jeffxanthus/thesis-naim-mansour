function e = end(o,k,n)
% tomSym/end - Overloaded - Get the last index of a tomSym object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if n==1
    e = o.s(end).sz1*o.s(end).sz2;
else
    if k==1
        e = o.s(end).sz1;
    else
        e = o.s(end).sz2;
    end
end

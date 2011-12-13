function y=rem(a,b)
% tomSym/rem - Overload the remainder function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if isequal(a,b)
    y = spzeros(size(a));
else
    y = quickop(mfilename,a,b);
end

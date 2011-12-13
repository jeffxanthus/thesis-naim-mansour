function y=plus(a,b)
% tomArrayIdx/plus - Add an offset to an index

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if isnumeric(a) && numel(a)==1
    y = b;
    y.i = a + b.i;
elseif isnumeric(b) && numel(b)==1
    y = a;
    y.i = a.i + b;
else
    error('Wrong type offset for tomArrayIdx');
end

function y = operand(n,o)
% tomSym/operand - Get an operand from a tomSym
%
% y = operand(n,o) gets the nth operand from tomSym object o.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if n==1 && length(o.s(end).a)==1 && o.s(end).a>0
    y = o;
    y.s(end) = [];
    y.sh = [];
    y.ih = [];
else
    ni = o.s(end).a(n);
    y = subsymb(ni,o);
end

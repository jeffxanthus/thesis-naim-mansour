function d = hownonlinear(f)
% hownonlinear - overloaded function
%
% d = hownonlinear(f) returns:
% d = 0, if f is a constant
% d = 1, if f is linear in all variables that it contains
% d = 2, if f is quadratic
% d > 2 if none of the above is true.
%
% Sometimes hownonlinear will return a higher value that is the actual
% case, because it failed to completely analyze the expression.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-08 by rutquist for TOMLAB release 7.7

sp = zeros(1,length(f.s));

for i=1:length(sp)
    subp = zeros(1,length(f.s(i).a));
    for k=1:length(subp)
        ix = f.s(i).a(k);
        if ix>0;
            subp(k) = sp(ix);
        else
            subp(k) = 0;
        end
    end
    sp(i) = getNonlinearity(i,f,subp);
end

d = sp(end);

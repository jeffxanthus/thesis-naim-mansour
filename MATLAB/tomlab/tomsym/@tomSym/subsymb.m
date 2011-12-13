function y = subsymb(ni,o)
% tomSym/subsymb - Get a sub-symbol from a tomSym
%
% y = subsymb(ni,o) gets the nth sub-symbol from tomSym object o.
%
% If ni<0, then a constant is returned.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-10-25 by rutquist for TOMLAB release 7.7

if(ni>0)
    % Extract a subset of the symbols list
    y = o;
    offs = length(y.d)+1;
    used = false(1,length(y.s)+offs);
    used(ni+offs) = true;
    for i = ni:-1:1
        if used(i+offs)
            used(y.s(i).a+offs) = true;
        end
    end
    ix = find(used)-offs;
    ixp = ix(ix>0);
    ixm = fliplr(-ix(ix<0));
    y.s = y.s(ixp);
    y.d = y.d(ixm);
    remap = zeros(size(used));
    remap(-ixm+offs) = -1:-1:-length(ixm);
    remap(ixp+offs) = 1:length(ixp);
    for i=1:length(y.s)
        y.s(i).a = remap(y.s(i).a+offs);
    end
    y.sh = [];
    y.ih = [];
else
    y = o.d{-ni};
end

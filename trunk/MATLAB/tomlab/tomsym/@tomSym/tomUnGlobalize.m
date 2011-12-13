function o = tomUnGlobalize(o)
% tomSym/tomUnGlobalize - undo the effect of tomGlobalize
%
% o = tomUnGlobalize(o)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

for i=1:length(o.s)
    if strcmp(o.s(i).op,'tomUnGlobalize')
        o.s(i).op = 'constant';
        o.d{-o.s(i).a(1)} = tomUnGlobalize(o.d{-o.s(i).a(1)},o.d{-o.s(i).a(2)},o.d{-o.s(i).a(3)});
        o.s(i).a = o.s(i).a(1);
    end
end

function o = squeeze(o)
% tomArray/squeeze
% Remove singelton dimensions

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-26 by rutquist for TOMLAB release 7.7

nonsing = (o.sz~=1);
o.sz = o.sz(nonsing);
o.ni = o.ni(nonsing(1:length(o.ni)));

function y = srepmat(a,sz1,sz2)
% tomSym/srepmat - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% $Id$

if numel(sz1)~=1 || numel(sz2)~=1
    error('Size argument to repmat must be scalar');
end

if numel(a)~=1
    y = repmat(a,sz1,sz2);
end

% Simplify if size = [ 1 1]
if sz1==0 || sz2==0
    y = zeros(sz1,sz2);
elseif sz1==1 && sz2==1
    y = a;
else
    y = tomSym(mfilename, sz1, sz2, a, sz1, sz2);
end


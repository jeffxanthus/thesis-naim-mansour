function y = maxIndicator(a,dim)
% maxIndicator - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin~=2
    dim = 1;
end

if size(a,dim)==1
    % Max of exactly one element
    y = ones(size(a));
else
    y = tomSym(mfilename,size(a,1),size(a,2),a,dim);
end

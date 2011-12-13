function [m,n] = size(o,i)
% tomCmplx/size - Get the size of a tomSym object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

if nargout>1
    [n,m] = size(o.re);
else
    if nargin==2
        m = size(o.re,i);
    else
        m = size(o.re);
    end
end

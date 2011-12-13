function o = reshape(o,siz,varargin)
% tomArray/reshape - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin>2
    siz = [siz horzcat(varargin{:})];
end

if prod(o.sz)~=prod(siz)
    error('Number of elements must not change.');
end

o.sz = siz;
o.ni = {};

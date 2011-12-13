function y = subs(y,varargin)
% tomArray/subs - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y.X = subs(y.X,varargin{:});

if isnumeric(y.X)
    y = unArray(y);
end


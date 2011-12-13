function z = zero_g(x,Prob) %#ok
% zero_g - The gradient of zero_f
%
% z = zero_g(x,Prob) returns an all-zero vector with numel(x) elements.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7


z = spalloc(1,numel(x),0);

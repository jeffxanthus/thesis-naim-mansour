function v = vec(M)
% vec - Transform a tomArray into a column vector
%
%  vec(M) is the same as M(:)
%
% For tomArrays, the result of vec is no longer a tomArray.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

v = vec(M.X);

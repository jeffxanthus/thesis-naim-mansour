function M = unName(M)
% unName - remove index names from a tomArray
%
%  X = unName(M) removes the index names from a tomArray, making it
%  possible to work with the array using conventional Matlab indexing, or
%  to assign new names.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

M.ni = {};

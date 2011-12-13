function p = pattern(f)
% pattern - The sparsity pattern of a matrix
%
% p = pattern(f) computes the sparse logical matrix p, that is
% true (1) wherever f is nonzero, and false (0) elsewhere.
%
% This function is overloaded for symbolic tomSym objects.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

p = sparse(f~=0);

function Y = scalerows(sc, M)
% Scalerows - Scale the rows of a matrix
%
% Y = scalerows(SC, M) is equivalent to Y = setdiag(SC)*M
%
% Each row in M is multiplied by the corresponding element in SC.
%
% This function is defined so that it can be overloaded for tomSym,
% reducing the occurence of setdiag whose derivative is very large and very
% sparse.
%
% See also submatrix

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

Y = setdiag(sc)*M;

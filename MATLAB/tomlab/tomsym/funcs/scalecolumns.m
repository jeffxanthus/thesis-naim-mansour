function Y = scalecolumns(sc, M)
% Scalecolumns - Scale the columns of a matrix
%
% Y = scalecolumns(SC, M) is equivalent to Y = M*setdiag(SC)
%
% Each column in M is multiplied by the corresponding element in SC.
%
% See also: scalerows, submatrix

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

Y = M*setdiag(sc);

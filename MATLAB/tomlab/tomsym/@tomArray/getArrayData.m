function X = getArrayData(M)
% unArray - extract data from a tomArray
%
%  X = getArrayData(M) extracts the data from tomArray M.
%
% Unlike unArray, getArrayData does not preserve the shape of M. The
% extracted data X will have the same number of elements as M, but it may
% have any shape.
%
% See also: tomArray/unArray tomArray/vec

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

X = M.X;

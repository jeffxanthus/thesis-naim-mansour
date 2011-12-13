function y = srepmat(a,sz1,sz2)
% tomSym/srepmat - implicit repmat (for tomSym internal use)
%
% srepmat(a,sz1,sz2) is the same as repmat(a,sz1,sz2) but tomSym uses this
% function to denote an implicit repmat, resulting from expressions like
% y+A, where y is a scalar and A is a matrix.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% $Id$

y = repmat(a,sz1,sz2);

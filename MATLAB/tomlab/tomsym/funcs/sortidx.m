function i = sortidx(x,dim,mode)
% sortidx - Second output from sort (indexes)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

[scratch,i] = sort(x,dim,mode);

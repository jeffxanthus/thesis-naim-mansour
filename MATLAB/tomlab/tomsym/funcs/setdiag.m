function D = setdiag(v)
% setdiag - Create a diagonal matrix from a column vector
%
% D = setidag(v) creates a sparse diagnoal n-by-n matrix D from the vector
% v with lenght n. If v is a matrix then it is first converted into a
% vector taking the elements columns-first.
%
% While setdiag() does approximately the same thing as diag(), it is more 
% clear what it does. (Since diag() also does the job of getdiag() when 
% presented with a matrix.)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-11-20 by rutquist for TOMLAB release 7.7

D = diag(sparse(v(:)));


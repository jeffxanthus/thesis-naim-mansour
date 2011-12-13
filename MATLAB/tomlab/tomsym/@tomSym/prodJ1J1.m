function J = prodJ1J1(M,dim)
% tomSym/prodJ1J1 - overloaded function 

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

J = tomSym(mfilename,size(M,3-dim)*numel(M),numel(M),M,dim);

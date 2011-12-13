function z = spzeros(m,n)
% spzeros - A sparse version of zeros.
%
% z = spzeros(n) gives an n-by-n matrix of zeros.
%
% z = spzeros(m,n) or spzeros([m,n]) is an m-by-n matrix of zeros.
%
% This is equivalent to sparse(m,n)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-01-24 by rutquist for TOMLAB release 7.7

if(nargin<2)
    if length(m)==2
        n = m(2);
        m = m(1);
    elseif length(m)==1
        n = m;
    else
        error('Illegal size for spzeros.');
    end
end

z = spalloc(m,n,0);

function s = prod(M,dim)
% tomSym/prod - Overloaded function 

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    if size(M,1)~=1
        dim = 1;
    else
        dim = 2;
    end
end

if ~isnumeric(dim) && (dim==1 || dim==2)
    error('Dimension argument must be either 1 or 2');
end

if dim==1
    sz = [1, size(M,2)];
else
    sz = [size(M,1), 1];
end

if prod(sz)==1 && tomCmp(M,'ctranspose')
    % Switch "dim" argument, rather than doing a transpose.
    dim = 3-dim;
    M = M';
end

if size(M,dim)==0
    % Empty product
    s = ones(sz);
elseif size(M,dim)==1
    % Product of exactly one element
    s = M;
elseif size(M,dim)==2
    % Product of two elements
    if dim==1
        s = submatrix(M,1,':').*submatrix(M,2,':');
    else % dim == 2
        s = submatrix(M,':',1).*submatrix(M,':',2);
    end
else
    s = tomSym(mfilename,sz(1),sz(2),M,dim);
end

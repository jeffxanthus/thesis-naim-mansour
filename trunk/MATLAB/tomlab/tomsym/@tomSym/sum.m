function s = sum(M,dim)
% tomSym/sum - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-22 by rutquist for TOMLAB release 7.7

if nargin<2
    if isempty(M)
        s = 0;
        return
    end
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
    % Empty sum
    s = spzeros(sz);
elseif size(M,dim)==1
    % Sum of exactly one element
    s = M;
elseif isdiag(M)
    if dim==1
        s = getdiag(M)';
    else
        s = getdiag(M);
    end
elseif tomCmp(M,'scalerows')
    if dim==1
        s = vec(operand(1,M))'*operand(2,M);
    else
        s = scalerows(operand(1,M),sum(operand(2,M),dim));
    end
elseif tomCmp(M,'scalecolumns')
    if dim==2
        s = operand(2,M)*vec(operand(1,M));
    else
        s = scalecolumns(operand(1,M),sum(operand(2,M),dim));
    end
else
        
    s = tomSym(mfilename,sz(1),sz(2),M,dim);
end

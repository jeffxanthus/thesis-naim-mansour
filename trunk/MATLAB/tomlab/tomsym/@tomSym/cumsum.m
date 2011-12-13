function y = cumsum(a,dim)
% tomSym/cumsum - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if isempty(a)
    y = a;
    return
end

if nargin<2
    if size(a,1)==1 && size(a,2)>1
        dim = 2;
    else
        dim = 1;
    end
end

if dim==1
    y = tril(ones(size(a,1)))*a;
else % dim==2
    y = a*triu(ones(size(a,2)));
end

function y = diff(x,n,dim)
% tomSym/diff - Difference - Overloaded function
%
% DIFF(X,N,DIM) givers the Nth difference of X along dimension DIM.
% This is the tomSym version of the built-in function with the same name.
%
% Diff does not compute the derivative! - See tomSym/derivative for that.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2 || isempty(n);
    n = 1;
end
if nargin<3
    if size(x,1)>1
        dim = 1;
    else
        dim = 2;
    end
end

if ~isnumeric(n) || ~isnumeric(dim)
    % TODO - This should be fixed when we alow variable size symbolic matrices.
    if tomCmp(n,'tom')
        disp('Did you mean to call "derivative"?');
    end
    error('The N and DIM arguments must be constants.');
end

nn = size(x,dim);
if nn<=n;
    % More differences than elements - return empty matrix.
    switch dim
        case 1
            y = zeros(0,size(x,2));
        case 2
            y = zeros(size(x,1),0);
        otherwise
            error('DIM must be equal to 1 or 2.');
    end
else
    % Return answer as a multiplication by a constant rather than a "diff"
    % symbol.
    A = speye(nn);
    for i=1:n
        A = spdiags([-ones(nn-i,1) ones(nn-i,1)],[0 1],nn-i,nn-i+1)*A;
    end
    switch dim
        case 1
            y = A*x;
        case 2
            y = x*A';
        otherwise
            error('DIM must be equal to 1 or 2.');
    end
end

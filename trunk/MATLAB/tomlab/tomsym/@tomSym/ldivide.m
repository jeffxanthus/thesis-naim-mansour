function y=ldivide(a,b)
% tomSym/ldivide - Overload the .\ left division operator
%
% The .\ operator can mean either matrix-matrix, scalar-matrix, or
% matrix-scalar division

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(a)==1
    % Matrix divided by scalar - easily translated into a product
    y = inv(a)*b;
elseif all(size(a) == size(b))
    % Matrix-matrix element-wise division
    sz = size(a);

    % We don't check for division-by zero. It will trigger a warning anyway.
    % Simplify if one of the arguments is zero or one.
    if iszero(b)
        y = zeros(sz);
    elseif isallone(a)
        y = b;
    elseif isequal(a,b)
        y = ones(sz);
    else
        y = quickop(mfilename,a,b);
    end
elseif numel(b)==1
    % Matrix divided into scalar
    y = (a.^-1)*b;
else
    error('Matrix dimensions must agree.');
end

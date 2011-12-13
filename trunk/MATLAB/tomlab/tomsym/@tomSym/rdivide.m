function y=rdivide(a,b)
% tomSym/rdivide - Overload the ./ division operator
%
% The ./ operator can mean either matrix-matrix, scalar-matrix, or
% matrix-scalar division

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if numel(b)==1
    % Matrix divided by scalar - easily translated into a product
    y = inv(b)*a;
elseif all(size(a) == size(b))
    % Matrix-matrix element-wise division
    sz = size(a);

    % We don't check for division-by zero. It will trigger a warning anyway.
    % Simplify if one of the arguments is zero or one.
    if iszero(a)
        y = zeros(sz);
    elseif isallone(b)
        y = a;
    elseif isequal(a,b)
        y = ones(sz);
    else
        y = quickop(mfilename,a,b);
    end
elseif numel(a)==1
    % Scalar divided by matrix
    y = a*b.^-1;
else
    error('Matrix dimensions must agree.');
end


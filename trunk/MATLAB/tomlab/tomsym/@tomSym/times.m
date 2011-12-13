function y=times(a,b)
% tomSym/times - Overload the .* multiplication operator (Hadamard product)
%
% The .* operator can mean either matrix-matrix, scalar-matrix, or
% matrix-scalar multiplication. We map them as follows:
%   times   - matrix-matrix element-wise multiplication
%   smtimes - scalar-matrix multiplication, where a is a scalar. It is
%             assumed that this operator commutes, so matrix-scalar
%             multiplication is mapped to this, with the order of operands
%             changed.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-26 by rutquist for TOMLAB release 7.7

% Move out minus signs
if tomCmp(a,'uminus')
    y = -(-a.*b);
    return
end
if tomCmp(b,'uminus')
    y = -(a.*-b);
    return
end

if all(size(a) == size(b))
    % Matrix-matrix element-wise multiplication
    sz = size(a);
    if numel(a)==1
        % Map scalar multiplication to mtimes instead
        y = a*b;
    elseif iszero(a) || iszero(b)
        y = zeros(sz);
    elseif isallone(a)
        y = b;
    elseif isallone(-a)
        y = -b;
    elseif isallone(b)
        y = a;
    elseif isallone(-b)
        y = -a;
    elseif isrepmat(a) && numel(b)>1
        y = smtimes(lookup(a,1),b);
    elseif isrepmat(b) && numel(a)>1
        y = smtimes(lookup(b,1),a);
    else
        y = quickop(mfilename,a,b);
    end
elseif numel(a)==1
    % Scalar-matrix multiplication
    y = smtimes(a,b);
elseif numel(b)==1
    % Scalar-matrix multiplication
    y = smtimes(b,a);
else
    error('Matrix dimensions must agree.');
end

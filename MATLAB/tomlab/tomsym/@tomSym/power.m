function y=power(a,b)
% tomSym/power - Overload the mpower operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if numel(b)==1
    % Matrix raised to a sclar, element-wise.
    y = spower(a,b);
elseif all(size(a) == size(b))
    % Element-wise power
    if iszero(a) && iszero(b)
        error('0^0 not allowed with tomSym');
    elseif iszero(a)
        y = zeros(sz);
    elseif iszero(b)
        y = ones(sz);
    elseif isallone(a)
        y = a; % 1^x = 1
    elseif isallone(b)
        y = a; % x^1 = x
    elseif tomCmp(a,'power') && isnumeric(b) && alll(b==round(b))
        y = power(operand(1,a),operand(2,a).*b);
    elseif tomCmp(a,'abs') && isnumeric(b) && alll(b/2==round(b/2));
        y = power(operand(1,a),b);
    else
        y = quickop(mfilename,a,b);
    end
elseif numel(a)==1
    y = repmat(a,size(b)).^b;
else
    error('Matrix dimensions must agree.');
end

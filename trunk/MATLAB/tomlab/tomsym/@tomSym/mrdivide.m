function y=mrdivide(a,b)
% tomSym/mrdivide - Overload the divide operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Error if the denominator is zero
if(iszero(b))
    error('tomSym divide by zero');
end

if numel(b)==1
    y = a./b;
else
    if size(a,2) ~= size(b,2)
        error('Matrix dimensions must agree.');
    end

    sz1 = size(a,1);
    sz2 = size(b,1);

    if iszero(a)
        % Multiplication by zero
        y = zeros(sz1,sz2);
    elseif isone(b)
        % Division by identity matrix
        y = a;
    elseif isnumeric(b) && isone(-b)
        % Division by negative identity matrix
        y = -a;
    elseif tomCmp(b,'setdiag')
        % Diagonal matrix is inverted element-wise
        y = a*inv(b);
    else
        y = tomSym(mfilename,sz1,sz2,a,b);
    end
end

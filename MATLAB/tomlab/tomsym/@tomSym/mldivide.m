function y=mldivide(a,b)
% tomSym/mrdivide - Overload the left divide operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Error if the denominator is zero
if(iszero(a))
    error('tomSym divide by zero');
end

if numel(a)==1
    y = a.\b;

else
    if size(a,1) ~= size(b,1)
        error('Matrix dimensions must agree.');
    end

    sz1 = size(a,2);
    sz2 = size(b,2);

    if iszero(b)
        % Multiplication by zero
        y = zeros(sz1,sz2);
    elseif isone(a)
        % Division by identity matrix
        y = b;
    elseif isnumeric(a) && isone(-a)
        % Division by negative identity matrix
        y = -b;
    elseif tomCmp(b,'setdiag')
        % Diagonal matrix is inverted element-wise
        y = inv(a)*b;
    else
        y = tomSym(mfilename,sz1,sz2,a,b);
    end
end

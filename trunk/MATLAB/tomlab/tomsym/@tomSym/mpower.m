function y=mpower(a,b)
% tomSym/mpower - Overload the mpower operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if size(a,1)~=size(a,2)
    error('Matrix power is only defined for square matrices.');
end
if numel(b)~=1
    error('TomSym Matrix power is only defined for scalar exponent.');
end
if numel(a)==1
    y = spower(a,b);
elseif iszero(a) && iszero(b)
    error('0^0 not allowed with tomSym');
elseif iszero(a)
    y = zeros(size(a));
elseif iszero(b)
    y = speye(size(a));
elseif isone(a)
    y = a; % 1^x = 1
elseif isone(b)
    y = a; % x^1 = x
elseif tomCmp(a,'mpower') && isnumeric(b) && b==round(b)
    y = mpower(operand(1,a),operand(2,a)*b);
elseif isnumeric(b) && b<0 && b==round(b)
    y = inv(a)^-b;
else
    y = tomSym(mfilename,size(a,1),size(a,2),a,b);
end

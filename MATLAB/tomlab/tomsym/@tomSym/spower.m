function y=spower(a,b)
% tomSym/spower - Element-wise power to a scalar exponent.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-18 by rutquist for TOMLAB release 7.7

if numel(b)~=1
    error('Spower expects a scalar exponent.');
end

if isempty(a)
    y = a;
elseif isnumeric(a) && nnz(a)<numel(a)
    if iszero(b)
        error('0^0 not allowed with tomSym');
    end
    [i,j,v] = find(a);
    y = sparse(i,j,spower(v,b),size(a,1),size(a,2));
elseif tomCmp(a,'setdiag')
    y = setdiag(spower(operand(1,a),b));
elseif tomCmp(a,'sparse')
    y = sparse(operand(1,a),operand(2,a),spower(operand(3,a),b),operand(4,a),operand(5,a));
elseif iszero(b)
    y = ones(size(a));
elseif isone(b)
    y = a; % x^1 = x
elseif tomCmp(a,'spower') && isnumeric(b) && b==round(b)
    y = spower(operand(1,a),operand(2,a)*b);
elseif tomCmp(a,'abs') && isnumeric(b) && alll(b/2==round(b/2));
    y = power(operand(1,a),b);
elseif isallone(a)
    y = a; % 1^x = 1
else
    y = tomSym(mfilename,size(a,1),size(a,2),a,b);
end

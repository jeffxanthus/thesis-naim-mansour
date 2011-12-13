function y=minus(a,b)
% tomSym/minus - Overload the minus operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-18 by rutquist for TOMLAB release 7.7

if numel(a)==1 && numel(b)~=1
    y = smplus(a,-b);
elseif numel(b)==1 && numel(a)~=1
    y = smplus(-b,a);
elseif tomCmp(a,'repmat') && numel(operand(1,a))==1
    y = operand(1,a)-b;
elseif tomCmp(b,'repmat') && numel(operand(1,b))==1
    y = a-operand(1,b);
elseif iszero(b)
    y = a;
elseif iszero(a)
    y = -b;
elseif isnumeric(a) && numel(a)>1 && all(size(a)==size(b)) && nnz(diff(a(:)))==0
    y = smplus(a(1),-b);
elseif isnumeric(b) && numel(b)>1 && all(size(a)==size(b)) && nnz(diff(b(:)))==0
    y = smplus(-b(1),a);
elseif isequal(a,b);
    y = spzeros(size(a));
elseif tomCmp(b,'uminus')
    y = a+operand(1,b);
elseif tomCmp(a,'plus') && isequal(operand(1,a),b)
    y = operand(2,a);
elseif tomCmp(a,'plus') && isequal(operand(2,a),b)
    y = operand(1,a);
    % TODO: There should be a more general varargin-sum.
else
    y = quickop(mfilename,a,b);
end

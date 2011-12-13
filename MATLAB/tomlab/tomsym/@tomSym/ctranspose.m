function p = ctranspose(a)
% tomSym/ctranspose - Overloaded operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-03 by rutquist for TOMLAB release 7.7

if numel(a)==1
    p = a;
elseif tomCmp(a, 'ctranspose')
    p = operand(1,a);
elseif tomCmp(a, 'setdiag') || tomCmp(a, 'setSymmetric')
    p = a;
elseif tomCmp(a, 'sparse')
    p = sparse(operand(2,a),operand(1,a),operand(3,a),operand(5,a),operand(4,a));
elseif tomCmp(a, 'repmat')
	p = repmat(operand(1,a)',operand(3,a),operand(2,a));
else
    p = tomSym(mfilename,size(a,2),size(a,1),a);
end

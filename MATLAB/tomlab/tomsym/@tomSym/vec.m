function y = vec(a)
% tomSym/vec - Transform a tomSym object into a column vector

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-02-08 by rutquist for TOMLAB release 7.7

% Simplify if a is already a column vector.
if size(a,2)==1
    y = a;
elseif tomCmp(a,'reshape')
    y = vec(operand(1,a));
elseif tomCmp(a,'repmat') && size(operand(1,a),2)==1
    y = repmat(operand(1,a),operand(2,a)*operand(3,a),1);
elseif tomCmp(a,'sparse')
    i = sub2ind(size(a),operand(1,a),operand(2,a));
    y = sparse(i,ones(size(i)),operand(3,a),numel(a),1);
elseif tomCmp(a,'ctranspose') && size(a,1)==1
    y = operand(1,a);
else
    y = tomSym('vec', numel(a), 1, a);
end

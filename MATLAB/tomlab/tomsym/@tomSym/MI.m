function y = MI(a)
% tomSym/MI - Matrix Inequality
%
% MI(A>=B) where A and/or B is a tomSym object is shorthand for
% positiveSemidefinite(A-B)
%
% MI(A>=b) where b is a scalar is shorthand for MI(A>=b*eye(size(A)),
% that is: All eigenvalues of A should be greater than or equal to b.
%
% Outside of MI(), A>=B would be an element-wise comparison.
%
% See also: positiveSemidefinite

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-29 by rutquist for TOMLAB release 7.7

if ~(tomCmp(a,'le') || tomCmp(a,'ge'))
    error('Matrix inequality must involve either >= or <=.');
end

if size(operand(1,a),1)~=size(operand(1,a),2)
    error('Matrix must be square.');
end

if numel(a) == 1
    y = a;
    return
end

o1 = operand(1,a);
o2 = operand(2,a);
if tomCmp(o1,'srepmat');
    o1 = setdiag(repmat(operand(1,o1),size(o1,1),1));
end
if tomCmp(o2,'srepmat');
    o2 = setdiag(repmat(operand(1,o2),size(o2,1),1));
end

if strcmp(operator(a),'le')
    y = positiveSemidefinite(-o1+o2);
else
    y = positiveSemidefinite( o1-o2);
end

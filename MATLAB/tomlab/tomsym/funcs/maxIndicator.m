function y = maxIndicator(a,dim)
% maxIndicator - Indicate the largest element
%
% maxIndicator(a,dim) returns a matrix of the same size as a, where the
% elements that would be selected by max(a,dim) are one, and the rest are
% zero.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin~=2
    dim = 1;
end

[ignore, i] = max(a,[],dim); %#ok
if dim==1
    y = sparse(i,1:size(a,2),ones(1,size(a,2)),size(a,1),size(a,2));
else
    y = sparse(1:size(a,1),i,ones(1,size(a,1)),size(a,1),size(a,2));    
end

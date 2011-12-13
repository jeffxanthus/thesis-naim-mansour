function y = subkron(a,b,ix1,ix2)
% subkron - submatrix of a Kronecker product.
%
% y = subkron(a,b,ix1,ix2) is equivalent to c = kron(a,b); y = c(ix1,ix2)
%
% Using subkron instead of kron can save a lot of memory in cases where
% only a small part of the matrix resulting from the Kronecker product is
% actually used.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if ischar(ix1) && strcmp(ix1,':')
    ix1 = 1:size(a,1)*size(b,1);
end
if ischar(ix2) && strcmp(ix2,':')
    ix2 = 1:size(a,2)*size(b,2);
end

ix1 = vec(ix1)';
ix2 = vec(ix2)';

ixa1 = vec(repmat(1:size(a,1),size(b,1),1))';
ixa2 = vec(repmat(1:size(a,2),size(b,2),1))';
ixb1 = vec(repmat((1:size(b,1))',1,size(a,1)))';
ixb2 = vec(repmat((1:size(b,2))',1,size(a,2)))';

ixa1 = ixa1(ix1);
ixa2 = ixa2(ix2);
ixb1 = ixb1(ix1);
ixb2 = ixb2(ix2);

y = a(ixa1,ixa2).*b(ixb1,ixb2);

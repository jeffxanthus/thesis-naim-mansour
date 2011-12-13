function y=ekron(M,n,p,A,B)
% tomSym/ekron - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% $Id$

if nargin < 4
    A = 'I';
end
if nargin < 5
    B = 'I';
end

if ~isnumeric(p) || (p~=1 && p~=2)
    error('p must equal 1 or 2');
end

if isone(A)
    A = 'I';
end
if isone(B)
    B = 'I';
end

if strcmp(A,'I')
    sz1 = size(M,1)*n;
    A0 = speye(sz1);
else
    sz1 = size(A,1);
    A0 = A;
end

if strcmp(B,'I')
    sz2 = size(M,2)*n;
    B0 = speye(sz2);
else
    sz2 = size(B,2);
    B0 = B;
end

% TODO: Replace A0 and B0 by if/else statements below.

if iszero(M) || iszero(A) || iszero(B)
    y = zeros(sz1,sz2);
elseif numel(M)==1
        y = A0*(M*speye(n))*B0;
elseif n==1
    y = A0*M*B0;
elseif size(A0,1) == 1
    if p==1
        y = vec(M'*reshape(A,size(M,1),n))'*B0;
    else
        y = vec(reshape(A,n,size(M,1))*M)'*B0;
    end
elseif size(B0,2) == 1
    if p==1
        y = A0*vec(M*reshape(B,size(M,2),n));
    else
        y = A0*vec(reshape(B,n,size(M,2))*M');
    end
else
    y = tomSym(mfilename,sz1,sz2,M,n,p,A,B);
end

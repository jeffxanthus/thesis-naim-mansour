function y=ekron(M,n,p,A,B)
% ekron - Kronecker product with an identity matrix.
%
% y = ekron(M,n,p,A,B) is equivalent to:
%
%   if p==1
%      y = A*kron(eye(n),M)*B
%   elseif p==2
%      y = A*kron(M,eye(n))*B
%   end
%
% The matrices A and B can be replaced by 'I' to represent an identity 
% matrix. For example: ekron(M,n,1,'I','I') is the same as kron(eye(n),M).
%
% If A is a row vector, or if B is a column vector, then the Kronecker
% product is never computed explicitly.
%
% TomSym uses ekron for reasons of memory efficiency.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% $Id$

if nargin < 4
    A = 'I';
end
if nargin < 5
    B = 'I';
end

if p==1
    if strcmp(A,'I')
        if strcmp(B,'I')
            y = kron(speye(n),M);
        else
            if size(B,2) == 1;
               y = vec(M*reshape(B,size(M,2),n));  
            else
                y = kron(speye(n),M)*B;
            end
        end
    elseif strcmp(B,'I')
        if size(A,1) == 1
            y = vec(M'*reshape(A,size(M,1),n))';
        else
            y = A*kron(speye(n),M);
        end
    else
        if size(B,2) == 1
            y = A*vec(M*reshape(B,size(M,2),n));
        elseif size(A,1) == 1
            y = vec(M'*reshape(A,size(M,1),n))'*B;
        else
            y = A*kron(speye(n),M)*B;
        end
    end
elseif p==2
    if strcmp(A,'I')
        if strcmp(B,'I')
            y = kron(M,speye(n));
        else
            y = sparse(size(M,1)*n,size(B,2));
            M = M';
            for i=1:size(B,2)
                y(:,i) = vec(reshape(B(:,i),n,size(M,1))*M);
            end
        end
    elseif strcmp(B,'I')
        if size(A,1) == 1
            y = vec(reshape(A,n,size(M,1))*M)';
        else
            y = A*kron(M,speye(n));
        end
    else
        if size(A,1) == 1
            y = vec(reshape(A,n,size(M,1))*M)'*B;
        else
            y = sparse(size(A,1),size(B,2));
            M = M';
            for i=1:size(B,2)
                y(:,i) = A*vec(reshape(B(:,i),n,size(M,1))*M);
            end
        end
    end
else
    error('p must equal 1 or 2.');
end

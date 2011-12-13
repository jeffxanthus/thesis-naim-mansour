% Back substitution of an upper-triangular matrix
%
%      solve A * x = b
%
% INPUT:
% A is n x n upper-triangular matrix
% b is n x 1 vector
%
% OUTPUT
% x is n x 1 vector
%
% function  x = backsubR(A,b)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 2, 1995.   Last modified Oct 12, 1999.

function  x = backsubR(A,b)

n    = length(b);
if n == 0
   x=[];
   return
end

detA = A(n,n);
x    = zeros(n,1);
x(n) = b(n)/A(n,n);
for j = n-1:-1:1,
    detA = detA*A(j,j);
    if detA == 0, break, end
    x(j) = (b(j) - A(j,j+1:n)*x(j+1:n))/A(j,j);
end
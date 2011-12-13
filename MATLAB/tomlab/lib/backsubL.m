% Back substitution of a lower-triangular matrix
%
%      solve A * x = b
%
% INPUT:
% A is n x n lower-triangular matrix
% b is n x 1 vector
%
% OUTPUT
% x is n x 1 vector
%
% function  x = backsubL(A,b)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Aug 8, 1999.      Last modified Oct 12, 1999.

function  x = backsubL(A,b)

n    = length(b);

if n==0
   x = [];
   return
end

detA = A(1,1);
x    = zeros(n,1);
x(1) = b(1)/A(1,1);
for j = 2:n,
    detA = detA*A(j,j);
    if detA == 0, break, end
    x(j) = (b(j) - A(j,1:j-1)*x(1:j-1))/A(j,j);
end
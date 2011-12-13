% AppRowQR:
%
% Append row to QR decomposition
%
% Given    A = Q * R, add row x to A which gives [A;x(:)']
%
% AppRowQR is using Givens rotations to compute the updated QR factorization
%
% function [R, Q] = AppRowQR(x, R, Q) 
%
%
% INPUT PARAMETERS
% x        Additional vector to be added as the new row in A
% Q,R      Matrices in QR factorization  A = Q * R 
% 
% OUTPUT PARAMETERS
%
% Q,R      Matrices in the expanded QR factorization  [A;x(:)'] = Q * R 
%
% If Q is not wanted, AppRowQR avoids the work of computing Q.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Sep 10, 1999. Last modified Jan 17, 2005.

function [R, Q] = AppRowQR(x, R, Q) 

% Append row last in QR decomposition

if nargin < 3
   error('AppRowQR needs three parameters x, R and Q !!!');
end

Full=nargout > 1;

n = size(Q,1);
p = size(R,2);

R=[R;x(:)'];

if Full
   Q = [Q zeros(n,1);zeros(1,n) 1];
end

I = [1;n+1];

for k=1:p
    I(1)=k;
    G=givens(R(k,k),R(n+1,k));
    R(I,k:p)=G*R(I,k:p);
    if Full
       Q(:,I)=Q(:,I)*G';
    end
end

% MODIFICATION LOG
%
% 050117 med   mlint revision
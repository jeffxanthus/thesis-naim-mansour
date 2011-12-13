% Compute QR factorization of matrix A and the pseudo rank pRank
% determined from the rank tolerance epsRank
%
% function [Q, R, E, pRank] = ComputeQR(A, epsRank) 
%
% INPUT PARAMETERS
% A        The matrix A
% epsRank  Rank tolerance. Default 1E-8
% 
% OUTPUT PARAMETERS
%
% Q,R,E    Matrices in QR factorization  A = Q * R * E'
% pRank    Pseudo rank

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written Aug 24, 1999. Last modified Aug 13, 2009.

function [Q, R, E, pRank] = ComputeQR(A, epsRank) 

if nargin < 2
   epsRank=1E-8;
end

[m n]= size(A);

if m > 1000
   ASPARSE=1;
else
   ASPARSE=0;
end

if ASPARSE
   P = colamd(A);
   [Q R]   = qr(A(:,P));  
   E = sparse(P,(1:n),ones(n,1),n,n);
   if n > 1 && m > 1
      maxR  = max(abs(diag(R)));
      pRank = full(nnz(abs(diag(R)) >= maxR*epsRank));
   else
      maxR  = abs(R(1,1));
      pRank = 1;
   end
else
   [Q R E] = qr(full(A));
   if ~isempty(R), maxR=abs(R(1,1)); 
   else 
       maxR=0; 
   end
   if n > 1 && m > 1
      pRank = nnz(abs(diag(R)) >= maxR*epsRank);
   else
      pRank = 1;
   end
   E=sparse(E);
end
if maxR == 0
   pRank=0;
end

% MODIFICATION LOG:
%
% 050726  med  colmmd call changed to colamd
% 090813  med  mlint check
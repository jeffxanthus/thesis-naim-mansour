% InsertQR:
%
% Insert column, variable VarNumber, into QR factorization of matrix A,
% at position pRank + 1
%
% function [Q, R, P, pRank, maxR, absR] = InsertQR(Q, R, P, pRank, ...
%           VarNumber, VarColumn, epsRank) 
%
% INPUT PARAMETERS
% Q,R       Matrices in QR factorization  A = Q * R * E'
% P         P is permutation made by E, [P, jE] = find(sparse(E));
% pRank     Pseudo rank
% VarNumber The variable number for the new variable
% VarColumn The matrix column for the new variable VarNumber
% epsRank   Rank tolerance
% 
% OUTPUT PARAMETERS
%
% Q,R      Matrices in QR factorization  A = Q * R * E'
% P        P is permutation made by E, [P, jE] = find(sparse(E));
% pRank    Pseudo rank
% maxR     Maximal diagonal element on diagonal of R.
% absR     abs(diag(R))

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Aug 24, 1999. Last modified Sept 26, 2000.

function [Q, R, P, pRank, maxR, absR] = InsertQR(Q, R, P, pRank, ...
          VarNumber, VarColumn, epsRank) 

% Insert column ?, variable VarNumber, into QR-decomposition, at pos pRank+1
%
% Maximum m = size(Q,1) columns are returned in P.

nargin;

[Q,R]=qrinsert(Q,R,pRank+1,VarColumn);

%for i=length(P):-1:size(R,2)
%    % Drop number i. Reorder P
%    j = P(i); 
%    P(P > j)=P(P > j) - 1;
%end

if ~isempty(P)
   %if length(VarNumber)~=1
   %   VarNumber
   %end
   ix=P >= VarNumber;
   if ~isempty(ix)
      P(ix)=P(ix) + 1;
   end
   %P(P >= VarNumber)=P(P >= VarNumber) + 1;
end


P=[P(1:pRank);VarNumber;P(pRank+1:length(P))];

if any(size(R) == 1)
   absR = abs(R(1,1));
else
   absR = abs(diag(R));
end

maxR  = max(absR);
pRank = nnz(absR >= maxR*epsRank);

%if size(R,2) > size(R,1)
%   % Remove extra unnecessary columns
%   for i=size(R,2):-1:size(R,1)+1
%       % Drop number i. Renumber P
%       j = P(i); 
%       P(P > j)=P(P > j) - 1;
%   end
%   P=P(1:size(R,1));
%   R=R(:,1:size(R,1));
%end
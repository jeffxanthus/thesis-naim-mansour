% DeleteQR:
%
% Delete the column in the QR decomposition corresponding to variable
% VarNumber.
%
% function [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, VarNumber)
%
% INPUT PARAMETERS
% Q,R       Matrices in QR factorization  A = Q * R * E'
% P         P is permutation made by E, [P, jE] = find(sparse(E));
% pRank     Pseudo rank
% VarNumber The variable number for the variable to delete
%
% OUTPUT PARAMETERS
%
% Q,R      Matrices in QR factorization  A = Q * R * E'
% P        P is permutation made by E, [P, jE] = find(sparse(E));
% pRank    Pseudo rank

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Aug 24, 1999. Last modified Aug 24, 1999.

function [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, VarNumber)

j=find(VarNumber==P);
if j <= pRank
    pRank=pRank-1;
end
[Q,R]=qrdelete(Q,full(R),j);
P=[P(1:j-1);P(j+1:length(P))];
P(P > VarNumber)=P(P > VarNumber) - 1;

% MODIFICATION LOG:
%
% 040728 med  Removed pragma
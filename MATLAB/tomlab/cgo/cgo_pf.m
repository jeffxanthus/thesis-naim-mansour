% function Fp = cgo_pf(F, c, Prob)
%
% Compute penalty function based on costly function f(x) 
% (and costly constraint mC-vector Cc(x) if defined)
% for set of n column vectors F and mC x n matrix c
%
% INPUT 
% f     n costly function values f(x_i), i=1,...,n
% c     mC x n matrix with costly constraint values 
%       c(j, x_i), j=1,...,mC, i=1,...,n
%
% OUTPUT
% Fp    Exact penalty function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 1, 2008.   Last modified July 1, 2008.

function Fp = cgo_pf(F, c, Prob)

if Prob.simType > 0
   n  = length(F);
   ixEQ = find(Prob.c_L == Prob.c_U);
   ixL  = find(~isinf(Prob.c_L) & Prob.c_L ~= Prob.c_U);
   ixU  = find(~isinf(Prob.c_U) & Prob.c_L ~= Prob.c_U);
   nEQ  = length(ixEQ);
   nL   = length(ixL);
   nU   = length(ixU);
   P    = zeros(n,1);
   if nEQ > 0 
      P = P + sum(abs(c(ixEQ,:)),1)';
   end
   if nL > 0 
      P = P + sum(max(0,Prob.c_L(ixL)*ones(1,n)-c(ixL,:)),1)';
   end
   if nU > 0
      P = P + sum(max(0,c(ixU,:)-Prob.c_U(ixU)*ones(1,n)),1)';
   end
   Fp = F + Prob.CGO.PEN*P;
else
   Fp = F;
end

% MODIFICATION LOG
%
% 080701 hkh  Written

% Sep_rJ.m
%
% Implements Algorithm II in Ruhe-Wedin (1980) for
% Separable Nonlinear Least Squares problem
% The problem is minimize f(x,z), where x is outer variables and
% z is inner (hidden) variables.
%
% In each iteration x^k the search step in x is computed as:
%       Js p = - rS, where rS and Js are corrected using the solution of
% min f(x^k,z) over z
%
% Let Jz=df(x^k,z)/dz at the minimal z value.
%
% Then [Q R] = qr(Jz). Let p be the rank of Jz and m the number of residuals,
% then the residuals r and Jacobian J are projected to the orthogonal
% complement of Jz, i.e.
%
%     rS = Q(:,p+1:m)' * r
%     Js = Q(:,p+1:m)' * J
%
% If Prob.LS.SepAlg == 1 this function is called
%
% function [rS,Js,Q,R,pRank] = Sep_rJ(x, r, J, z, Jz, Prob)
%
% x     Current iterate
% r     Current weighted residual
% J     The Jacobian matrix wrt x
% Jz    The Jacobian matrix wrt z

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written Oct 25, 1998.  Last modified Aug 13, 2009.

function [rS,Js,Q,R,pRank] = Sep_rJ(x, r, J, z, Jz, Prob)

global LS_A

DEBUG=0;

epsRank=Prob.optParam.eps_Rank;

[Q R]=qr(full(Jz));

m=size(Jz,1);
n=size(J,2);

if any(size(R) == 1)
   absR  = abs(R(1,1));
   pRank = absR > 0;
else
   absR=abs(diag(R));
   pRank=sum(absR > epsRank*abs(R(1,1))); 
end

if pRank~=size(R,2) & DEBUG
   fprintf('Rank of Q in qr(Jz) %d\n',pRank);
end

% Compute projection onto orthogonal complement to R(Jz)
rS=[Q(:,pRank+1:m)'*r;zeros(pRank,1)];

if DEBUG
   fprintf('Weighted residual r\n');
   xprint(r,'r:  ');
   fprintf('Projected residual Q^T * rW\n');
   xprint(rS,'rS: ');
   xprinte(r-rS,'rS-r');
end

% Compute the same projection for the Jacobian matrix
Js=[Q(:,pRank+1:m)'*J;zeros(pRank,n)];

% Check what is left (the error in the projection:

rSerr=Q(:,1:pRank)'*r;

if sum(rSerr) >=1E-12 & DEBUG
   fprintf('Sum of what is left, should be close to 0: %20.17f\n',sum(rSerr));
   fprintf('Or near the accuracy of the data\n');
   xprinte(Q(:,1:pRank)'*r);
end

% MODIFICATION LOG:
%
% 981025  hkh  Routine written using old code in exp_r.m
% 981108  hkh  Add input parameter z, the solution in the separated subproblem
% 000926  hkh  Safeguard diag command for any size == 1
% 090813  med  mlint check
% exp_J.m
%
% function J = exp_J(x, Prob);
%
% Computes Jacobian to least squares problem formulation of
% the fitting of positive sums of exponentials to empirical data 
%
% If Prob.LS.SepAlg==1
% use separable nonlinear least squares formulation,
% computing alpha (and beta) as the solution to a NNLS (non-negative linear
% least squares) problem
% This computation is done in exp_r.m and alpha,beta (eType==4) is sent as 
% extra parameters in Prob structure, Prob.alpha and Prob.beta
% to this routine. If not alpha, beta sent, they are computed here as well.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written May 21, 1995.  Last modified Jun 5, 2008.

function J = exp_J(x, Prob)

global LS_A wLS

y  = Prob.LS.y;
t  = Prob.LS.t;
m  = size(y,1);

[SepAlg, p, wType, eType] = expGet1(Prob);

nE = p + p*(eType==5);

if SepAlg
   E     = LS_A;
   if isfield(Prob,'alpha')
      alpha = Prob.alpha;
      beta  = Prob.beta;
   else
      % Call is not from exp_r, compute alpha (and beta)
      [alpha,beta,Jz,wLS]=expLS(x,E,Prob);
   end
   % Only J=dr/d(lambda) is computed, lambda = b. Lambda is hidden.
   switch eType
    case 1
      J=-(t*alpha').*E;
    case 2
      J=(t*alpha').*E;
    case 3
      J=(t.^2*(-alpha)').*E;
    case 4
      %beta=x(nE+p+1:nE+2*p);
      J=-(t*beta'+t.^2*alpha').*E;
    case 5
      % Also compute J=dr/d(my)
      lam=x(1:p);
      my=x(p+1:2*p);
      J=[(t*ones(1,p)-ones(m,1)*my').*E*diag(alpha), E*diag(-alpha.*lam)];
    otherwise
   end
else
   alpha=x(nE+1:nE+p);
   E=LS_A;
   switch eType
    case 1
      J=[-(t*alpha').*E,     E];
    case 2
      J=[(t*alpha').*E, 1-E];
    case 3
      J=[(t.^2*(-alpha)').*E,  diag(t)*E];
    case 4
      %ix =[nE+p+1:nE+2*p];
      %beta=x(ix);
      beta=x(nE+p+1:nE+2*p);
      J=[-(t*beta'+t.^2*alpha').*E,  diag(t)*E, E];
    case 5
      lam=x(1:p);
      my=x(p+1:2*p);
      J=[(t*ones(1,p)-ones(m,1)*my').*E*diag(alpha), E*diag(-alpha.*lam),1-E];
    otherwise
   end
end

% MODIFICATION LOG:
%
% 981018 hkh Get Yt from Prob.NLLS.Yt, not Prob.Yt. Same for t.
% 981022 hkh Changes for new NLPLIB, global variables revised.
% 981025 hkh Rewrite this routine
% 981102 hkh Error in eType==2 corrected
% 990216 hkh Expand routine to handle eType=5; distribution estimation
% 011205 hkh Use Prob for alpha, beta. Use switch statement
% 011205 hkh Safeguard, compute alpha & beta, when exp_J not called from exp_r
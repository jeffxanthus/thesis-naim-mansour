% exp_r.m
%
% function [r, J, z, Jz] = exp_r(x, Prob)
%
% Computes residuals to Nonlinear Least Squares problem formulation of
% the fitting of positive sums of exponentials to empirical data
%
% OUTPUT:
% r   residual
% J   Jacobian
%
%     If Prob.LS.SepAlg == 1 use separable non linear least squares
%     formulation, computing alpha (and beta if eType==4) as the solution to
%     a NNLS problem.  (NNLS = non-negative linear least squares problem)
% z   alpha (and beta if eType==4)
% Jz  Derivative w.r.t. (alpha,beta)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 21, 1995.  Last modified Dec 5, 2001.

function [r, J, z, Jz] = exp_r(x, Prob)

global LS_A wLS

% Note: b = lambda, a = alpha, c = beta, d = my
[SepAlg, p, wType, eType] = expGet1(Prob);

y  = Prob.LS.y;
t   = Prob.LS.t(:);
m   = size(y,1);
lambda=max(0,x(1:p)); % Take max

if eType <= 4
   E=exp(t*(-lambda)');
   nE=p;
elseif eType == 5
   % 2 dimensional exp function
   nE=p+p;
   my=x(p+1:nE);
   E = exp((t*ones(1,p)-ones(m,1)*my')*diag(-lambda));
end

LS_A=E;

if SepAlg
   % Determine alpha by NNLS. Jz is weighted by wLS
   [alpha,beta,Jz,wLS]=expLS(x,E,Prob);
   z=[alpha;beta];
   Prob.alpha=alpha;
   Prob.beta=alpha;
   % Send alpha and beta in Prob to Jacobian computing.
   % wLS is sent as global
   %J=nlp_J([x;z],Prob);
   J=nlp_J(x,Prob);
else
   z = []; J = []; Jz= [];
   alpha=x(nE+1:nE+p);
   if eType == 4
      beta=x(nE+p+1:nE+2*p);
   end
end

switch eType
  case 1
   % eType = 1, sum(j=1:m) [ sum(i=1:p) a(i) exp(-b(i)t(j)) - y(j) ]^2
   yMod=E*alpha;   
  case 2
   % eType = 2, sum(j=1:m) [ sum(i=1:p) a(i) (1-exp(-b(i)t(j))) - y(j) ]^2
   yMod=(1-E)*alpha;   
  case 3
   % eType = 3, sum(j=1:m) [ sum(i=1:p) t(j) a(i) exp(-b(i)t(j)) - y(j) ]^2
   yMod=t.*(E*alpha);   
  case 4
   % eType = 4, sum(j=1:m) [ sum(i=1:p) (t(j)a(i)-c(i)) exp(-b(i)t(j))-y(j)]^2
   yMod=t.*(E*alpha)+E*beta;
  case 5
   % Distribution estimation
   % eType = 5, sum(j=1:m) [ sum(i=1:p) (a(i) exp(-b(i)(t(j)-d(i))-y(j)]^2
   yMod=(1-E)*alpha;
end
r=yMod-y;

% It is dangerous to avoid adding y, as the weighting and separation fails
% if not the same transformation are made on y.
 
% MODIFICATION LOG:
%
% 981018  hkh  Get y from Prob.NLLS.y, not Prob.y. Same for t.
% 981022  hkh  Changes for new NLPLIB, global variables revised.
% 981024  hkh  Totally revised the routine. New logic in LS handling.
% 981025  hkh  Fix handling of Separable NLLS
% 981102  hkh  Error in comments for eType==2
% 981108  hkh  Add variable z as output
% 981118  hkh  NOT use new flag Prob.NLLS.yUse, to determine if to add Yt
% 990216  hkh  Expand routine to handle eType=5; distribution estimation
%              Will not work for separable algorithms
% 011205  hkh  Use Prob for alpha, beta. Use switch statement
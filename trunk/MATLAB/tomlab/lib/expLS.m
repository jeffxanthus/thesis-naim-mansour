% expLS.m
%
% Computes linear least squares solution to
% the fitting of positive sums of exponentials to empirical data 
% when the exponential parameters has fixed values.
% Used in separable formulation
%
% function [alpha, beta, Jz, wLS] = expLS(x, E, Prob)
%
% INPUT
% x      x parameters
% E      Exponential matrix
% Prob   TOMLAB problem structure
%
% OUTPUT
% alpha  Linear weights
% beta   Additional linear weights
% Jz     The jacobian wrt alpha (or alpha/beta)
% wLS    Weighting vector/matrix

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Oct 26, 1995.   Last modified Jan 16, 2003.

function [alpha, beta, Jz, wLS] = expLS(x, E, Prob)

[SepAlg, p, wType, eType] = expGet1(Prob);

y  = Prob.LS.y;
t   = Prob.LS.t;

% Determine alpha by non-negative least squares
r=[];
wLS=ExpFitW(x,r,Prob);
if eType==1
   Jz=diag(wLS)*E;
elseif eType==2;
   Jz=diag(wLS)*(1-E);
elseif eType==3;
   Jz=diag(wLS.*t)*E;
elseif eType==4;
   Jz=[diag(wLS.*t)*E,diag(wLS)*E];
elseif eType==5;
   Jz=diag(wLS)*(1-E);
end

%alpha=nnls1(Jz,wLS.*y);
alpha=Tnnls(Jz,wLS.*y);

if eType==4
   beta = alpha(p+1:p+p);
   alpha= alpha(1:p);
else
   beta = [];
end

% MODIFICATION LOG:
%
% 981026  hkh  Extract this code from exp_J
% 981211  hkh  Return wLS as fourth parameter
% 990216  hkh  Expand routine to handle eType=5; distribution estimation
% 030116  hkh  Change name to Tnnls
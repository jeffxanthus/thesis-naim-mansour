% expSetLa sets alpha and lambda parameters for exponential fitting in
% problem structure Prob.ExpFit
%
% function Prob = expSetLa(Prob, lambda, alpha, beta);
%
% INPUT:
%  Prob     Problem structure
%           p is assumed defined in Prob.ExpFit.p
%  lambda   p-dimensional lambda vector (intensities)
%  alpha    p-dimensional alpha  vector (weights)
%  beta     p-dimensional beta   vector (2nd weights when eType==4);
%
% OUTPUT:
%  Prob     Problem structure, changed entries Prob.ExpFit.alpha and
%           Prob.ExpFit.lambda.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 19, 1998.  Last modified June 21, 1999.

function Prob = expSetLa(Prob, lambda, alpha, beta)

if nargin < 4, beta=[];  end
if nargin < 3, alpha=[];  end
if nargin < 2, lambda=[]; end

p=Prob.ExpFit.p;

if ~isempty(beta)
   if length(beta)==p
      Prob.ExpFit.beta=beta; 
   end
else
  Prob.ExpFit.beta=[];
end

if ~isempty(alpha)
   if length(alpha)==p
      Prob.ExpFit.alpha=alpha; 
   end
else
  Prob.ExpFit.alpha=[];
end

if ~isempty(lambda)
   Prob.ExpFit.lambda=lambda; 
else
  Prob.ExpFit.lambda=[];
end
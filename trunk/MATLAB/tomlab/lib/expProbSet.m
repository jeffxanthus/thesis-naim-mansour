% expProbSet is setting values in the structure Prob. 
% Used by problem definition file exp_prob.m
%
% function Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
%                          qType, sigType, lambda, alpha, beta, x0Type, sumType)
%
% OUTPUT:
% Prob   Structure is set to these values.
               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written June 24, 1999. Last modified Nov 5, 2000.

function Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
                         qType, sigType, lambda, alpha, beta, x0Type, sumType)

if checkType('exp',Prob.probType)
   Prob.ExpFit.p=p;
   Prob.ExpFit.wType=wType;
   Prob.ExpFit.eType=eType;
   Prob.ExpFit.infCR=infCR;
   Prob.ExpFit.dType=dType;
   Prob.ExpFit.geoType=geoType;
   Prob.ExpFit.qType=qType;
   Prob.ExpFit.sigType=sigType;
   Prob.ExpFit.lambda=lambda;
   Prob.ExpFit.alpha=alpha;
   Prob.ExpFit.beta=beta;
   Prob.ExpFit.x0Type=x0Type;
   Prob.ExpFit.sumType=sumType;

else
   Prob.ExpFit=[];
end

% MODIFICATION LOG:
%
% 990624  hkh  Making routine from ProbSet
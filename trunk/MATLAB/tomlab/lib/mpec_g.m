% function g = mpec_g(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of the gradient f(x) when
% solving a complementarity problem with extra slack variables.
%
% mpec_g calls the original user routine for the gradient with the
% original variables only. Zeros are added for the slack variables as
% these do not enter into the objective.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function g = mpec_g(x,Prob,varargin)

g = feval(Prob.orgProb.FUNCS.g,x(1:Prob.orgProb.N),Prob.orgProb,varargin{:});

% Slacks don't enter in objective, just add zeros. 
g(end+1:Prob.N)=0.0;

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814 med  FUNCS used for callbacks instead
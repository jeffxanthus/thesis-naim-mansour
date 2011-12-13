% function f = mpec_f(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of function values f(x) when
% solving a complementarity problem with extra slack variables.
%
% mpec_f calls the original user routine for the objective function with the
% original variables only.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function f = mpec_f(x,Prob,varargin)

f = feval(Prob.orgProb.FUNCS.f,x(1:Prob.orgProb.N),Prob.orgProb,varargin{:});

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814 med  FUNCS used for callbacks instead
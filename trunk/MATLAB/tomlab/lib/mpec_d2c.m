% function d2c = mpec_d2c(x, lam, Prob, varargin)
%
% TOMLAB gateway routine for computation of the constraint Hessian d2c(x)
% when solving a complementarity problem with extra slack variables.
%
% mpec_d2c calls the original user routine for the constraint Hessian with
% the original variables only. If slack variables have been added to the
% problem, extra zero rows and columns are added.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function d2c = mpec_d2c(x,lam,Prob,varargin)

d2c = feval(Prob.orgProb.FUNCS.d2c,x(1:Prob.orgProb.N),lam,Prob.orgProb,varargin{:});

% If there are slacks, expand the constraint Hessian with zeros
% because the slacks enter only linearly in the constraints.

if(Prob.MPEC.ns>0)
   d2c(Prob.N,Prob.N) = 0.0;
end

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814 med  FUNCS used for callbacks instead
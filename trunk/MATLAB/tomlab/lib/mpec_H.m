% function H = mpec_H(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of the Hessian f(x) when
% solving a complementarity problem with extra slack variables.
%
% mpec_H calls the original user routine for the Hessian with the
% original variables only. If slack variables have been added to the
% problem, extra zero rows and columns are added.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function H = mpec_H(x,Prob,varargin)

H = feval(Prob.orgProb.FUNCS.H, x(1:Prob.orgProb.N), Prob.orgProb , varargin{:} );

% Slacks don't enter in objective. Just set a zero at the last slack
% variable diagonal to expand the matrix, if needed.

if (Prob.MPEC.ns>0)
    H(Prob.N,Prob.N) = 0.0;
end

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814 med  FUNCS used for callbacks instead
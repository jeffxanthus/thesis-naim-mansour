% function c = mpec_c(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of nonlinear constraints c(x) when
% solving a complementarity problem with extra slack variables.
%
% mpec_c calls the original user routine for the nonlinear constraints 
% with the original variables only.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function cx = mpec_c(x,Prob,varargin)

cx = feval(Prob.orgProb.FUNCS.c,x(1:Prob.orgProb.N),Prob.orgProb,varargin{:});

x=x(:);

% Subtract the slack variable values
if ~isempty(Prob.MPEC.MP)
   s  = x(Prob.orgProb.N+1:Prob.N);
   cx = cx - Prob.MPEC.MP*s;
end

% Duplicated constraints - added at the end of cx
if ~isempty(Prob.MPEC.exC)
   cx = [ cx(:) ; cx(Prob.MPEC.exC) ];
end

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814 med  FUNCS used for callbacks instead
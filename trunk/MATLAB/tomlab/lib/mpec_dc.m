% function dc = mpec_dc(x, Prob, varargin)
%
% TOMLAB gateway routine for computation of the constraint Jacobian dc(x)
% solving a complementarity problem with extra slack variables.
%
% mpec_dc calls the original user routine for the Jacobian with the
% original variables only. If slack variables have been added to the
% problem, their Jacobian coefficients are appended to the output.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 19, 2006.   Last modified Nov 7, 2011.

function dc = mpec_dc(x,Prob,varargin)

% Number of variables in original problem
n = Prob.orgProb.N;

% Evaluate original dc(x) function
Func = Prob.orgProb.FUNCS.dc;
if ~isempty(Func)
   dc = feval(Func,x(1:n),Prob.orgProb,varargin{:});
elseif ~isempty(Prob.orgProb.ConsPattern) & Prob.CheckNaN > 0
   % CheckNaN allows to have the original problem's Jacobian
   % estimated numerically.
   dc = Prob.orgProb.ConsPattern;
   dc(find(dc)) = NaN;
end

% Add pattern of slack variables, -1's according to Prob.MPEC.MP
dc = [dc,-Prob.MPEC.MP];

% If constraints have been duplicated, add these rows but only the 
% original problem columns
k = Prob.MPEC.exC(:);
if ~isempty(k)
   dc = [ dc ; dc(k,1:n),sparse(length(k),size(Prob.MPEC.MP,2)) ];
end

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060814  med FUNCS used for callbacks instead
% 111107 ango Interact with CheckNaN flag

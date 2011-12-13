% TOMLAB gateway routine
%
% nlp_d2c computes the 2nd part of the Hessian to the Lagrangian function,
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)
%
%
% nlp_d2c calls the routine Prob.FUNCS.d2c either as
%           d2c=feval(Prob.FUNCS.d2c, x, lam) or
%           d2c=feval(Prob.FUNCS.d2c, x, lam, Prob, varargin{:})
% depending on the number of inputs in Prob.FUNCS.d2c
%
% Prob.FUNCS.d2c returns lam' * d2c(x), i.e. if there are m constraints,
% d2c is the weighted sum of m matrices of size n by n. Each weight is
% the Lagrange parameter lam(i). This means that
%
% d2c = sum_i=1:m  lam(i) * d2c(i) / dx^2
%
% d2c(i) / dx^2 is the Hessian matrix (second order matrix) for constraint i
%
% function d2c=nlp_d2c(x, lam, Prob, varargin)
%
% The global counter variable n_d2c is incremented

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 10, 1998.   Last modified July 25, 2011.

function d2c=nlp_d2c(x, lam, Prob, varargin)

global n_d2c NLP_xdc NARG
global mad_dc
n_d2c=n_d2c+1; 

Func = Prob.FUNCS.d2c;
x=x(:);

N = min(length(x),Prob.N);

if Prob.ADCons == -1
   if length(NLP_xdc) ~= length(x)
       NLP_xdc = [];
   end
   if ~isempty(NLP_xdc) & all(NLP_xdc == x)
       z=getderivs(mad_dc);
       if isempty(z)
          d2c = sparse(Prob.N,Prob.N);
       else
          d2c=lam'*z(:,:);
          d2c=sparse(reshape(d2c,Prob.N,Prob.N));
       end
   else
       nlp_dc(x(1:N), Prob, varargin);
       z=getderivs(mad_dc);
       if isempty(z)
          d2c = sparse(Prob.N,Prob.N);
       else
          d2c=lam'*z(:,:);
          d2c=sparse(reshape(d2c,Prob.N,Prob.N));
       end
   end
elseif isempty(Func) | Prob.ConsDiff ~= 0
   % Here we call a numerical difference routine
   % First check if unconstrained
   if isempty(Prob.FUNCS.c)
      d2c=[];
   elseif abs(Prob.ConsDiff) == 7
      d2c=FDcHess_par(x(1:N),lam,Prob,[],varargin{:});
   else
      d2c=FDcHess(x(1:N),lam,Prob,[],varargin{:});
   end
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(6);
   end
   if p > 3
      d2c = feval(Func, x(1:N), lam, Prob, varargin{:});
   elseif p==3
      d2c = feval(Func, x(1:N), lam, Prob);
   elseif p==2
      d2c = feval(Func, x(1:N), lam); 
   else
      d2c = [];
   end
   if Prob.CheckNaN ~= 0
      [iN,jN,d2cN] = find(isnan(d2c));
      if ~isempty(iN)    % There are elements set to NaN, to be estimated
         % Only set one row with the variables needed to be estimated
         Prob.ConsPattern = sparse(ones(length(iN),1),...
              jN,d2cN,length(lam),size(d2c,2));
         d2cN=FDcHess(x(1:N),lam,Prob,[],varargin{:});
         % Merge analytic and numerical dc
         d2c(isnan(d2c)) = d2cN(isnan(d2c));
      end
   end
end

% MODIFICATION LOG:
%
% 981028  hkh  Remove ctrl vector for scaling
% 981120  hkh  Error in comments
% 981126  hkh  Use xnargin as filter, to avoid bug in Matlab5.1
% 990909  hkh  Redesign, avoid computing structure information.
% 020409  hkh  Use global NARG instead of calling xnargin every time
% 030113  ango Add check for numerical Hessian, calls FDcHess
% 030114  ango Fixed dangerous bug in call to FDcHess
% 030114  hkh  Not need returns deleted
% 030127  hkh  Check NaN elements, estimate numerically; Wrong ConsDiff test
% 030228  ango Changed way to get size of ConsPattern for CheckNaN option
% 031201  hkh  Revising AD handling, new for MAD.
% 031204  hkh  lam' * AD-matrix , not AD-matrix * lam, is OK
% 031206  hkh  Use getderivs for MAD. Empty implies 0 matrix
% 040526  hkh  Use x(1:N) in all function calls, define N
% 060814  med  FUNCS used for callbacks instead
% 070713  frhe lam' * AD-matrix(:,:) is OK. Added checks on NLP_xdc
% 110726  hkh   Added FDcHess (parfor) call if abs(Prob.NumDiff)==7

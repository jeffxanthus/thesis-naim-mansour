% TOMLAB gateway routine 
%
% nlp_d2r computes the 2nd part of the Hessian for a 
% nonlinear least squares problem, i.e.
%
%    
%           sum(i=1:m) r(i) * d2r_i / dx_j dx_k
%
% in
%
%   f(x) =   0.5 * r' * r, where r is a m-vector of residuals
%   g(x) =         J' * r, where J is the m by n Jacobian
%
% nlp_d2r calls the routine Prob.FUNCS.d2r either as 
%           d2r=feval(Prob.FUNCS.d2r, x ) or
%           d2r=feval(Prob.FUNCS.d2r, x, r) or
%           d2r=feval(Prob.FUNCS.d2r, x, r, J) or
%           d2c=feval(Prob.FUNCS.d2r, x, r, J, Prob, varargin{:})
% depending on the number of inputs in Prob.FUNCS.d2r
%
% function d2r=nlp_d2r(x, Prob, r, J, varargin)
%
% The global counter variable n_d2r is incremented

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Oct 23, 1998.  Last modified Aug 14, 2006.

function d2r=nlp_d2r(x, Prob, r, J, varargin)

global n_d2r LS_xJ
global mad_J

n_d2r=n_d2r+1;

if nargin < 4
  J = [];
  if nargin < 3
    r = [];
  end
end

if isempty(r)
  r = nlp_r(x, Prob);
end

if isempty(J)
  J = nlp_J(x, Prob);
end

x=x(:);
N = min(length(x),Prob.N);

if Prob.ADObj == -1
   % Hopefully multiplying by r from the right is correct, have to be tested
   if all(LS_xJ == x)
      d2r=getinternalderivs(mad_J)*r;
   else
      nlp_J(x(1:N), Prob, varargin);
      d2r=getinternalderivs(mad_J)*r;
   end
elseif isempty(Prob.FUNCS.d2r)
   if isempty(Prob.FUNCS.J)
      d2r=[];
   else   % Numerical differences possible
      d2r=FDrHess(x(1:N), Prob, r, J, varargin{:});
   end
   return
else
   n_d2r=n_d2r+1; 
   p=xnargin(Prob.FUNCS.d2r);
   if p>4
      d2r=feval(Prob.FUNCS.d2r, x(1:N), r, J, Prob, varargin{:});
   elseif p==4
      d2r=feval(Prob.FUNCS.d2r, x(1:N), r, J, Prob);
   elseif p==3
      d2r=feval(Prob.FUNCS.d2r, x(1:N), r, J); 
   elseif p==2
      d2r=feval(Prob.FUNCS.d2r, x(1:N), r); 
   elseif p==1
      d2r=feval(Prob.FUNCS.d2r, x(1:N)); 
   else
      d2r=[];
   end
end

% MODIFICATION LOG
%
% 981023 hkh  Routine written, based on nlp_d2c
% 981018 hkh  Delete ctrl vector for scaling
% 981126 hkh  Use xnargin as filter, to avoid bug in Matlab5.1
% 990909 hkh  Revised avoiding compuing structural information.
% 031201 hkh  Revising AD handling, new for MAD
% 040526 hkh  Use x(1:N) in all function calls, define N
% 050121 frhe Added num. diff. and input arg. order changed
% 060814 med  FUNCS used for callbacks instead
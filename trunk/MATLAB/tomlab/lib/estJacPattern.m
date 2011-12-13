% function JacPattern = estJacPattern(Prob,Trials,varargin)
%
% estJacPattern estimates Prob.JacPattern, the sparsity pattern of the
% residual Jacobian for nonlinear least squares problems
%
% The routine generates Trials random initial points between lower and
% upper bound
%
% It first tries to call the analytic Jacobian, if available.
% If no analytic Jacobian, it tries TOMLAB /MAD (if installed)
% Otherwise it estimates a numerical Jacobian
%
% INPUT:
% Prob        The Tomlab problem structure
% Trials      Number of trials to increase the probability that no element
%             by chance is 0 in the point tried. Default 2
%
% OUTPUT:
% JacPattern  A 0-1 m by n-matrix, sparse, with the residual Jacobian
%             pattern, see the description of Prob.JacPattern

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Apr 14, 2004.  Last modified Aug 13, 2009.

function JacPattern = estJacPattern(Prob,Trials,varargin)
if nargin < 2
   Trials = [];
end
if isempty(Trials), Trials = 2; end

global NARG
Func = Prob.FUNCS.r;
if isempty(NARG)
   p = xnargin(Func);
else
   p = NARG(7);
end
N    = Prob.N;
BIG  = 1000;
x_0  = Prob.x_0;
if isempty(x_0)
   x_0 = zeros(N,1);
end
M   = length(Prob.LS.y);
if M==0
   r = nlp_r(x_0,Prob,varargin{:});
   M = length(r);
end
x_L  = Prob.x_L(:);
nL   = length(x_L);
if nL < N
   x_L(nL+1:N) = x_0(nL+1:N)-BIG;
   x_L = x_L(:);
end
x_U  = Prob.x_U(:);
nU   = length(x_U);
if nU < N
   x_U(nU+1:N) = x_0(nU+1:N)+BIG;
   x_U = x_U(:);
end
% Safe guard x_0 inside bounds
x_0 = min(x_U,max(x_L,x_0));
xD              = x_U - x_L;
xD(isinf(xD))   = 2*BIG;
x_L(isinf(x_L)) = x_0(isinf(x_L))-BIG;

Prob.x_L = x_L;
Prob.x_U = x_U;

J = [];

if Prob.NumDiff == 0 & ~isempty(Prob.FUNCS.J)
   % Try analytic J
   J  = sparse(M,N);
   for i = 1:Trials
       x   = x_L+rand(N,1).*xD;
       J_x = abs(nlp_J(x,Prob,varargin{:}));
       if isempty(J_x)
          J = [];
          break
       else
          J = J + J_x;
       end
   end
elseif checkMAD(0)
   J  = sparse(M,N);
   % Try MAD for J
   try
      for i = 1:Trials
          x    = x_L+rand(N,1).*xD;
          if p > 2
             mad_r=feval(Func, fmad(x(1:N),speye(N)),Prob,varargin{:});
          elseif p > 1
             mad_r=feval(Func, fmad(x(1:N),speye(N)), Prob);
          else
             mad_r=feval(Func, fmad(x(1:N),speye(N)));
          end
          %r=getvalue(mad_r);
          J = J + abs(getinternalderivs(mad_r));
      end
   catch
       % Try FD for J
       Prob.NumDiff = 1;
       for i = 1:Trials
           x = x_L+rand(N,1).*xD;
           J = J + abs(nlp_J(x,Prob,varargin{:}));
       end
   end
end
if isempty(J)
   J  = sparse(M,N);
   % Try FD for J
   Prob.NumDiff = 1;
   for i = 1:Trials
       x = x_L+rand(N,1).*xD;
       J = J + abs(nlp_J(x,Prob,varargin{:}));
   end
end

JacPattern = spones(J);
% Set JacPattern 0 for all fixed variables
JacPattern(:,find(xD==0)) = 0;

% MODIFICATION LOG
%
% 040414  hkh  Algorithm formulated and written, based on estConsPattern
% 040728  hkh  Emergency call to nlp_r must have input x_0, not x
% 040901  med  getvalue lower case
% 041123  hkh  Safe guard for empty nlp_J Jacobian
% 050901  med  Corrected Trials default to 2, updated help
% 050901  med  nL was used instead of nU for x_U check
% 060609  hkh  Safe guard x_0 inside [x_L,x_U], Pattern 0 for fixed variables
% 060814  med  FUNCS used for callbacks instead
% 090813  med  mlint check
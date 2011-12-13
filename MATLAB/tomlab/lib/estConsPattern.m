% function ConsPattern = estConsPattern(Prob,Trials,varargin)
%
% estConsPattern estimates Prob.ConsPattern, the sparsity pattern of the
% constraint Jacobian
%
% The routine generates Trials random initial points between lower and
% upper bound
%
% It first tries to call the analytic constraint Jacobian, if available.
% If no analytic Jacobian, it tries TOMLAB /MAD (if installed)
% Otherwise it estimates a numerical constraint Jacobian
%
% INPUT:
% Prob        The Tomlab problem structure
% Trials      Number of trials to increase the probability that no element
%             by chance is 0 in the point tried. Default 2
%
% OUTPUT:
% ConsPattern A 0-1 m by n-matrix, sparse, with the constraint
%             Jacobian pattern, see the description of Prob.ConsPattern

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Apr 13, 2004.  Last modified Aug 13, 2009.

function ConsPattern = estConsPattern(Prob,Trials,varargin)
if nargin < 2
   Trials = [];
end
if isempty(Trials), Trials = 2; end

N    = Prob.N;
M    = Prob.mNonLin;
if M == 0
   ConsPattern = [];
   return;
end

global NARG
Func = Prob.FUNCS.c;
if isempty(NARG)
   p = xnargin(Func);
else
   p = NARG(4);
end
dc   = [];
BIG  = 1000;
x_0  = Prob.x_0;
if isempty(x_0)
   x_0 = zeros(N,1);
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

if Prob.ConsDiff == 0 & ~isempty(Prob.FUNCS.dc)
   % Try analytic dc
   dc   = sparse(M,N);
   for i = 1:Trials
       x  = x_L+rand(N,1).*xD;
       dc_x=abs(nlp_dc(x,Prob,varargin{:}));
       if isempty(dc_x)
          dc = [];
          break;
       else
          dc = dc + dc_x;
       end
   end
elseif checkMAD(0)
   dc   = sparse(M,N);
   % Try MAD for dc
   try
      for i = 1:Trials
          x    = x_L+rand(N,1).*xD;
          if p > 2
             mad_c=feval(Func, fmad(x(1:N),speye(N)),Prob,varargin{:});
          elseif p > 1
             mad_c=feval(Func, fmad(x(1:N),speye(N)), Prob);
          else
             mad_c=feval(Func, fmad(x(1:N),speye(N)));
          end
          dc = dc + abs(getinternalderivs(mad_c));
      end
   catch
       % Try FD for dc
       Prob.ConsDiff = 1;
       for i = 1:Trials
           x  = x_L+rand(N,1).*xD;
           dc = dc + abs(nlp_dc(x,Prob,varargin{:}));
       end
   end
end

if isempty(dc)
   dc   = sparse(M,N);
   % Try FD for dc
   Prob.ConsDiff = 1;
   for i = 1:Trials
       x  = x_L+rand(N,1).*xD;
       dc = dc + abs(nlp_dc(x,Prob,varargin{:}));
   end
end

ConsPattern = spones(dc);
% Set ConsPattern 0 for all fixed variables
ConsPattern(:,find(xD==0)) = 0;

% MODIFICATION LOG
%
% 040413  hkh  Algorithm formulated and written
% 040901  med  getvalue lower case
% 040930  ango Return [] immediately if Prob.mNonLin is zero
% 041114  hkh  Safe guard for empty nlp_dc Constraint Jacobian
% 050901  med  Corrected Trials default to 2, updated help
% 050901  med  nL was used instead of nU for x_U check
% 060609  hkh  Safe guard x_0 inside [x_L,x_U], Pattern 0 for fixed variables
% 060814  med  FUNCS used for callbacks instead
% 081119  ango Check for no nonlinear constraints earlier
% 090813  med  mlint check
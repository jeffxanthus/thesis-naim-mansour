% function d2cPattern = estd2cPattern(Prob,Trials,varargin)
%
% estd2cPattern estimates Prob.d2cPattern, the sparsity pattern of the
% Hessian of the second part of the Lagrangian function
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)
%
% The routine generates Trials random initial points between lower and
% upper bound
%
% It first tries to call the analytic d2c routine, if available.
% If no analytic d2c, but analytic Jacobian routine,
% it tries TOMLAB /MAD (if installed)
% Otherwise it estimates a numerical finite difference d2c
%
% INPUT:
% Prob        The TOMLAB problem structure
% Trials      Number of trials to increase the probability that no element
%             by chance is 0 in the point tried. Default 2
%
% OUTPUT:
% d2cPattern  A 0-1 n by n-matrix, sparse, with the pattern of
%             d2c, see the description of Prob.d2cPattern

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sep 1, 2005.  Last modified Aug 13, 2009.

function d2cPattern = estd2cPattern(Prob,Trials,varargin)
if nargin < 2
   Trials = [];
end
if isempty(Trials), Trials = 2; end

global NARG
Funcdc  = Prob.FUNCS.dc;
Funcd2c = Prob.FUNCS.d2c;

if isempty(NARG)
   if isempty(Funcdc)
      pdc = 0;
   else
      pdc = xnargin(Funcdc);
   end
   if isempty(Funcd2c)
      pd2c = 0;
   else
      pd2c = xnargin(Funcd2c);
   end
else
   pdc  = NARG(5);
   pd2c = NARG(6);
end
N   = Prob.N;
BIG  = 1000;
x_0  = Prob.x_0(:);
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
   x_U(nU+1:N) = x_0(nL+1:N)+BIG;
   x_U = x_U(:);
end
% Safe guard x_0 inside bounds
x_0 = min(x_U,max(x_L,x_0));
xD              = x_U - x_L;
xD(isinf(xD))   = 2*BIG;
x_L(isinf(x_L)) = x_0(isinf(x_L))-BIG;

Prob.x_L = x_L;
Prob.x_U = x_U;

if Prob.NumDiff == 0 & pd2c > 0
   % Try analytic d2c
   d2c = sparse(N,N);
   for i = 1:Trials
       x = x_L+rand(N,1).*xD;
       d2c_x = abs(nlp_d2c(x,rand(Prob.mNonLin,1),Prob,varargin{:}));
       if isempty(d2c_x)
          d2c = [];
          break
       else
          d2c = d2c + d2c_x;
       end
   end
else
   d2c = [];
end
if isempty(d2c)
   d2c = sparse(N,N);
   if checkMAD(0)
      % Try MAD for d2c
      try
         for i = 1:Trials
             x    = x_L+rand(N,1).*xD;
             if pdc > 2
                mad_d2c=feval(Funcdc, fmad(x(1:N),speye(N)),Prob,varargin{:});
             elseif pdc > 1
                mad_d2c=feval(Funcdc, fmad(x(1:N),speye(N)),Prob);
             else
                mad_d2c=feval(Funcdc, fmad(x(1:N),speye(N)));
             end
             z=getderivs(mad_d2c);
             d2c = d2c + abs(sparse(reshape(rand(N,1)'*z,Prob.N,Prob.N)));
         end
      catch
          % Try FD for d2c
          if pdc > 0
             % Utilize analytic gradient routine
             Prob.ConsDiff = -1;
          else
             Prob.ConsDiff = 1;
          end
          for i = 1:Trials
              x = x_L+rand(N,1).*xD;
              nlp_c(x,Prob,varargin{:});
              nlp_dc(x,Prob,varargin{:});
              d2c = d2c + abs(nlp_d2c(x,rand(Prob.mNonLin,1),Prob,varargin{:}));
          end
      end
   else
      % Try FD for d2c
      if pdc > 0
         % Utilize analytic gradient routine
         Prob.ConsDiff = -1;
      else
         Prob.ConsDiff = 1;
      end
      for i = 1:Trials
          x = x_L+rand(N,1).*xD;
          nlp_c(x,Prob,varargin{:});
          nlp_dc(x,Prob,varargin{:});
          d2c = d2c + abs(nlp_d2c(x,rand(Prob.mNonLin,1),Prob,varargin{:}));
      end
   end
end

d2cPattern = spones(d2c);
% Set d2cPattern 0 for all fixed variables
d2cPattern(:,find(xD==0)) = 0;

% MODIFICATION LOG
%
% 050901  med  Written
% 050914  ango Fix a few incorrect nlp_d2c calls
% 060609  hkh  Safe guard x_0 inside [x_L,x_U], Pattern 0 for fixed variables
% 060814  med  FUNCS used for callbacks instead
% 090813  med  mlint check
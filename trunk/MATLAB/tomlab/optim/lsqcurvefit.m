% lsqcurvefit is the TOMLAB equivalent to LSQCURVEFIT in Optimization TB 
% The NLLS/NLP solver actually used is selectable, see Prob.Solver.Tomlab below.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% Opt tbx lsqcurvefit solves nonlinear least squares problems on the form:
%
%	min       0.5 * sum {r(x).^2} <==>      min       0.5 *  r(x)^T * r(x)
%    x                                       x
%             x_L <=   x  <= x_U                      x_L <=   x  <= x_U
%
% TOMLAB lsqcurvefit solves nonlinear least squares problems on the form:
%
%	min       0.5 * sum {r(x).^2} <==>      min       0.5 *  r(x)^T * r(x)
%    x                                       x
%             x_L <=   x  <= x_U                      x_L <=   x  <= x_U
%             b_L <=  Ax  <= b_U                      b_L <=  Ax  <= b_U
%
% where r is the m-dimensional residual vector defined as
%
%            r(x) = F(x,xData) - yData
%
% and x is a n-dimensional unknown parameter vector, with bounds x_L and x_U.
%
% The linear constraints are treated by the Tomlab solvers as standard
% but they are not treated in OPTIM TB 
% They could be added in the Tomlab Prob structure, input as extra
% argument or as global variable otxProb, see below
%
% function [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =...
%   lsqcurvefit(rFunc, x, xData, yData, x_L, x_U,  options, Prob, varargin)
%
% INPUT: ( 4 arguments always needed )
%
% rFunc    Function that computes the residual r(x) (and maybe J(x))
% x        Starting value for the design variables
% xData    Input  data in residual norm 0.5*||F(x,xData)-yData||^2
% yData    Output data in residual norm 0.5*||F(x,xData)-yData||^2
% x_L      Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% x_U      Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%          Jacobian, JacobPattern, LineSearchType, LevenbergMarquardt,
%          MaxFunEvals, MaxIter, DiffMinChange and DiffMaxChange,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%          JacobMult: If nonempty, Tomlab makes i=1:dim(x) calls
%             feval(JacobMult,J,e_i,-1,varargin),   to obtain the J matrix.
%          e_i= (....,1, ....), 1 in the ith position, 0:s otherwise
%          The Tomlab solvers cannot utilize the J*x computation in JacobMult
%
%    NOTE: If Jacobian is on, the user defines the Jacobian matrix
%          and returns it together with the residual vector in rFunc
%          [r_k,J_k]=feval(rFunc,x,Prob)
%          The Prob extra argument is possible in the TOMLAB lsqcurvefit,
%          and can be used to send information to rFunc.
%          If r_k is in R^m and x in R^n, then J_k is an m-by-n matrix
%          where J(i,j) is the partial derivative of r_k(i) with respect to x(j)
%
% Prob     The TOMLAB problem input structure, or the first extra user
%          input argument
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
% Additional fields used in the Prob structure, beside standard fields:
%
% Prob.Solver.Tomlab  Name of the Tomlab solver to use instead of the
%                     Optimization Toolbox solver LSQNONLIN
%
% Prob.Solver.Tomlab ='nlssol';   selects the Tomlab /SOL nlssol solver
% Prob.Solver.Tomlab ='clsSolve'; selects the Tomlab clsSolve solver
% Prob.Solver.Tomlab ='snopt';    selects the Tomlab /SOL snopt solver
% Prob.Solver.Alg 1-5 selects the 5 best different methods in clsSolve.
%
% Default is to use the solver returned by:
%            GetSolver('cls', length(x) > 200 | Prob.LargeScale,0)
%
% See the online help for cls§olve (help cls§olve) and nlssol (help nlssol and
% help nlssolTL. TL is added to name of interface file for all MEX solvers).
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.Solver.Tomlab = 'nlssol';
%
% OUTPUT:
%
% x        Optimal design parameters
% f_k      Optimal residual sum of squares, sum {r(x).^2} (Note! no 0.5)
% r_k      The optimal residual vector
% ExitFlag exit condition of lsqcurvefit.
%      > 0 lsqcurvefit converged to a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 Errors, ExitFlag=-Inform, see the Inform parameter description
%
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.funcCount     Number of function evaluations
%   Output.algorithm     Type of algorithm used
%   Output.steplength    Length of last step
%   Output.cgiterations  Number of CG iterations (if used)
%   Output.firstorderopt The first-order optimality (if used)
% Lambda   Structure with Lagrange multipliers at the solution
%    Lambda.lower     Lagrange multipliers for the lower bounds lb
%    Lambda.upper     Lagrange multipliers for the upper bounds ub
% J_k      The optimal Jacobian matrix
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format)
% Structure Result. Fields used (Also see help for clsSolve):
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   f_k      Function value 0.5*x'*F*x+c'*x
%   r_k      Residual at x_k
%   g_k      Gradient at x_k
%   J_k      Jacobian matrix at x_k
%   Solver   clsSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written July 28, 1999.   Last modified May 21, 2008.

function [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =...
  lsqcurvefit(rFunc, x, xData, yData, x_L, x_U,  options, Prob, varargin)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

if nargin < 8, Prob = [];
   if nargin < 7, options = [];
      if nargin < 6, x_U = []; 
         if nargin < 5, x_L = []; 
            if nargin < 4 
	       error('lsqcurvefit requires four input arguments rFunc and x');
end, end, end, end, end

% ------------Initialization----------------
global otxProb

if ~isempty(Prob) & ~isstruct(Prob)
   Psave=Prob;
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
   Prob.varargin = [{Psave},varargin];
elseif isempty(Prob)
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
   if nargin > 11
      Prob.varargin = [{[]},varargin];
   else
      Prob.varargin = varargin;
   end
else
   if isfield(Prob,'TOMLAB')
      Prob.varargin = varargin;
   else
      % Assume structure is part of extra input
      Psave=Prob;
      if ~isempty(otxProb)
         Prob = otxProb;
      else
         Prob = ProbDef;
      end
      Prob.varargin = [{Psave},varargin];
   end
end

Prob.LS.t    = xData(:);
Prob.LS.y    = yData(:);
if length(x) > 500 | length(Prob.LS.y) > 2000
   Prob.LargeScale = 1;
   if ~isfield(options,'LargeScale')
      options.LargeScale = 'on';
   elseif isempty(options.LargeScale)
      options.LargeScale = 'on';
   end
else
   Prob.LargeScale = 0;
   if ~isfield(options,'LargeScale')
      options.LargeScale = 'off';
   elseif isempty(options.LargeScale)
      options.LargeScale = 'off';
   end
end

% Avoid using Yt in TOMLAB NLLS test problems when defining the residual

Prob.LS.yUse = 0;

Prob.probType=checkType('ls');

Prob=ProbCheck(Prob,'lsqcurvefit',4);

if ~isempty(Prob.LS.t) & size(Prob.LS.t,1)~=length(Prob.LS.y)
   fprintf('Size xData: %d %d. Length yData %d\n',size(Prob.LS.t),length(Prob.LS.y));
   error('lsqcurvefit: Illegal sizes of xData (1st dim) and yData');
end

z = version;
M7 = z(1) == '7';
Prob.OPTTB.M7= 0;

if isa(rFunc,'cell')
   funArgIn  = 1;
   if length(rFunc) == 1
      Prob.FUNCSX.r=rFunc{1};
      funArgOut = -1;
   elseif length(rFunc) == 2
      Prob.FUNCSX.r=rFunc{1};
      Prob.FUNCSX.J=rFunc{2};
      if isempty(Prob.FUNCSX.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      else
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      end
   else
      Prob.FUNCSX.r=rFunc{1};
      Prob.FUNCSX.J=rFunc{2};
      Prob.FUNCSX.d2r=rFunc{3};
      if isempty(Prob.FUNCSX.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      elseif isempty(Prob.FUNCSX.d2r)
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      else
         funArgOut = -3;
         Prob.NumDiff=0;
      end
   end
elseif isa(rFunc,'function_handle')
   Prob.FUNCSX.r=rFunc;
   if M7
      Prob.OPTTB.M7 = 1;
      funArgIn  = abs(nargin(rFunc));
      funArgOut = abs(nargout(rFunc));
      if funArgOut == 1
         Prob.NumDiff = max(1,Prob.NumDiff);
      elseif funArgOut == 2
         Prob.NumDiff = min(-1,Prob.NumDiff);
      else
         Prob.NumDiff=0;
      end
   else
      SS = functions(rFunc);
      funArgIn  = abs(nargin(SS.function));
      funArgOut = abs(nargout(SS.function));
   end
else
   Prob.FUNCSX.r=rFunc;
   funArgOut = abs(nargout(rFunc));
   funArgIn  = abs(nargin(rFunc));
   if funArgOut == 1
      Prob.NumDiff = max(1,Prob.NumDiff);
   elseif funArgOut == 2
      Prob.NumDiff = min(-1,Prob.NumDiff);
   else
      Prob.NumDiff=0;
   end
end

Prob.FUNCSX.funArgOut=funArgOut;
Prob.FUNCSX.funArgIn =funArgIn;

Prob.FUNCSX.r=rFunc;
Prob.FUNCSX.x=x;
if funArgOut > 1
   options.Jacobian = 'on';
else
   options.Jacobian = 'off';
end

Prob.N = length(x(:));
Prob = rmfield(Prob,'TOMLAB');
Prob.varargin = [{Prob}, Prob.varargin];

[x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] =...
    lsqnonlin('ls2c_rJ', x, x_L, x_U,  options, Prob.varargin{:});

f_k=0.5*f_k;

Result.f_k=f_k;

Result.SolverAlgorithm = ['lsqcurvefit/ ' Result.SolverAlgorithm];

% MODIFICATION LOG
%
% 030129 hkh  Revision for v4.0. Use global variable otxProb
% 050422 hkh  Revision for Matlab 7. Set LargeScale.
% 060814 med  FUNCS used for callbacks instead
% 070526 med  Removed otxProb assignment
% 070708 hkh  Add check on options.LargeScale being []
% 080414 hkh  Revise help, add text about choice of solver
% 080521 med  Defaults call added
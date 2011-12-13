% fmincon is the TOMLAB equivalent to fmincon in Optimization TB
% The NLP solver actually used is selectable, see Prob.Solver.Tomlab below.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% TOMLAB fmincon solves generally constrained problems on the OPTTB form:
%
%	min                f(x)
%        x
%                 x_L <=   x        <= x_U
%                          A   * x  <= b
%                          Aeq * x  == beq
%                          c(x)     <= 0
%                          ceq(x)   == 0
%
% TOMLAB constrained solvers use the more general constrained problems format:
%
%	min              f(x)
%        x
%                 x_L <=   x   <= x_U
%                 b_L <=   Ax  <= b_U
%                 c_L <=  c(x) <= c_U
%
% where f(x) is the objective function and
% x is a n-dimensional unknown parameter vector, with bounds x_L and x_U.
% A is a mA by n matrix with linear constraints, with bounds b_L and b_U
% c(x) is a m vector of nonlinear constraints, with bounds c_L and c_U
%
% The OPTTB form is therefore transformed to
%
%	min              f(x)
%        x
%             x_L             <=   x   <= x_U
%             b_L = [-Inf ]   <=   A = [ A  ] * x  <= b_U = [ 0 ]
%                   [ beq ]   <=       [ Aeq]      <=       [beq]
%             c_L = [-Inf ]   <=   c = [ c(x)  ]   <= c_U = [ 0 ]
%                   [  0  ]   <=       [ ceq(x)]   <=       [ 0 ]
%
%
% function [x, f_k, ExitFlag, Output, Lambda, g, H, Result] = fmincon(...
%     Func, x_0, A, b, Aeq, beq, x_L, x_U, conFunc, options, Prob, varargin)
%
% INPUT: ( 2 arguments always needed )
%
% Func     Function that computes the function f(x) (maybe g(x),H(x))
% x_0      Starting value for the design variables
% A        Linear constraint matrix for inequalities
% b        Right hand side for linear inequality constraints
% Aeq      Linear constraint matrix for equalities
% beq      Right hand side for linear equality constraints
% x_L      Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% x_U      Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
% conFunc  Function that computes the constraints c(x) (maybe dc(x))
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%          GradObj, Hessian, HessPattern, HessUpdate, LineSearchType,
%          MaxFunEvals, MaxIter, DiffMinChange and DiffMaxChange,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%
%          The Prob extra argument is possible in the TOMLAB fmincon,
%          and can be used to send information to Func and conFunc.
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
%                     Optimization Toolbox solver FMINCON
%
% Prob.Solver.Tomlab ='snopt';    selects the Tomlab /SOL snopt solver
% Prob.Solver.Tomlab ='conSolve'; selects the Tomlab conSolve solver
%
% Default is to use the solver returned by:
%            GetSolver('con', length(x) > 200 | Prob.LargeScale,0)
%
% See the online help for conSolve (help conSolve) and snopt (help snopt and
% help snoptTL. TL is added to name of interface file for all MEX solvers).
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.Solver.Tomlab = 'snopt';
%
% Another important information for large-scale optimization is the
% 0-1 pattern of the Jacobian of the constraints.
% In the Tomlab format this matrix is defined as Prob.ConsPattern
% Set otxProb.ConsPattern, if sending a global Prob structure, to make sparse
% solvers like SNOPT utilize the pattern. It is very important for speed.
% Both in the case where analytic derivatives are given, and when numerically
% estimating them.
% Note that Prob.ConsPattern should be an m-by-n matrix, because Tomlab
% is using a different convention for how to define the Jacobian dc of the
% constraints:
% If c(x) in R^m and x in R^n, then dc is an m-by-n matrix where
% dc(i,j) is the partial derivative of c_i(x) with respect to x(j)
%
% OUTPUT:
%
% x        Optimal design parameters
% f_k      Final function value
% ExitFlag exit condition of fmincon.
%      > 0 fmincon converged to a solution X.
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
%
% Lambda   Structure with Lagrange multipliers at the solution
%    Lambda.lower        Lagrange multipliers for the lower bounds of x
%    Lambda.upper        Lagrange multipliers for the upper bounds of x
%    Lambda.ineqlin      Lagrange multipliers for the linear inequalities A
%    Lambda.eqlin        Lagrange multipliers for the linear equalities Aeq
%    Lambda.ineqnonlin   Lagrange multipliers for the nonlinear inequalities c(x)
%    Lambda.eqnonlin     Lagrange multipliers for the nonlinear equalities  ceq(x)
%
% g        The final gradient
% H        The final Hessian
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (Tomlab format)
% Structure Result. Fields used (Also see help for conSolve):
%   Iter     Number of iterations
%   ExitText Text about the Exit conditions
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   f_k      Function value at x_k
%   g_k      Gradient at x_k
%   H_k      Hessian matrix at x_k
%   Solver   Tomlab solver, e.g. SNOPT, KNITRO or conSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 28, 1999.  Last modified Aug 23, 2011.

function [x, f_k, ExitFlag, Output, Lambda, g, H, Result] = fmincon(...
    Func, x_0, A, b, Aeq, beq, x_L, x_U, conFunc, options, Prob, varargin)

global otxProb alphaV n_f

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

if nargin < 11, Prob = [];
   if nargin < 10, options = [];
      if nargin < 9, conFunc = [];
         if nargin < 8, x_U = []; 
            if nargin < 7, x_L = []; 
               if nargin < 6, beq = []; 
                  if nargin < 5, Aeq = []; 
                     if nargin < 4, b = []; 
                        if nargin < 3, A = []; 
            if nargin < 2 
	       error('fmincon requires two input arguments Func and x_0');
end, end, end, end, end, end, end, end, end, end

% ------------Initialization----------------

if ~isempty(Prob) & ~isstruct(Prob)
   Psave=Prob;
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
   %Prob.Psave=Psave;
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
z = version;
M7 = z(1) == '7';

solvType=checkType('con');
Prob.probType = DefPar(Prob,'probType',solvType);

n = max([length(x_0(:)),length(x_L(:)),length(x_U(:))]);
if n == 0
   error('x_0, x_L and x_U are empty');
end

if length(x_0) < n
    Prob.x_0 = zeros(n,1);
else
    Prob.x_0 = x_0(:);
end

if length(x_L) < n
    Prob.x_L = -inf*ones(n,1);
else
    Prob.x_L = x_L(:);
end

if length(x_U) < n
    Prob.x_U = inf*ones(n,1);
else
    Prob.x_U = x_U(:);
end

Prob.N       = n;
if isempty(x_0), x_0 = zeros(n,1); end

Prob.OPTTB.x     = x_0;
Prob.OPTTB.M7    = 0;
Prob.OPTTB.M7Con = 0;

mAeq=size(Aeq,1);
mA  =size(A,1);

b   = b(:);
beq = beq(:);

if mAeq==0
   Prob.A=A;
   Prob.b_L=-Inf*ones(mA,1);
   Prob.b_U=b;
elseif mA==0
   Prob.A=Aeq;
   Prob.b_L=beq;
   Prob.b_U=beq;
else
   Prob.A=[Aeq;A];
   Prob.b_L=[beq;-Inf*ones(mA,1)];
   Prob.b_U=[beq;b];
end
Prob.mLin = size(Prob.A,1);
Solver = DefPar(Prob.Solver,'Tomlab',[]);

if isempty(Solver)
   Solver = GetSolver('con',n > 200 | Prob.LargeScale,0);
   Prob.Solver.Tomlab = Solver;
end
% Hessian/gradient numerical approximation off
Prob.NumDiff=0; 
% Constraint Jacobian numerical approximation off
Prob.ConsDiff=0; 
% No Automatic Differentiation
Prob.ADObj=0;
Prob.ADCons=0;
% LargeScale on
Prob.LargeScale=1;
% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

Prob=ProbCheck(Prob,Solver,solvType);

% Set default options
% Diagnostics - off

Diagnostic=0;

% ??? Prob.optParam.DiffGradMaxChange=1E-1;

% Opt tbx defaults
% Display          final
% Diagnostics      off
% LargeScale       on
% MaxIter          400
% TolX             1E-6
% TolFun           1E-6
% TolCon           1E-6
% MaxFunEvals      100*numberOfVariables
% MaxIter          400
% LineSearchType   quadcubic, but in Tomlab Cubic is default
% DiffMinChange    1e-8

% Opt tbx defaults (not used)
% DiffMaxChange    1e-1 
% DerivativeCheck - off
% HessMult         []
% LevenbergMarq    on
% MaxIter          200
% PrecondBandWidth 0    
% TypicalX         ones(numberOfVariables,1) 
% TolPCG           0.1
% MaxPCGIter       max(1,floor(numberOfVariables/2))

if ~isempty(options)
   % Opt tbx defaults - test if changed
   optDef.TolFun      = 1E-6;
   optDef.TolX        = 1E-6;
   optDef.TolCon      = 1E-6;
   optDef.MaxIter     = 400;
   optDef.MaxFunEvals = 100*n; % MaxFunEvals 100*numberOfVariables
   % DiffMaxChange 1e-1 DiffMinChange 1e-8
   optDef.DiffMinChange = 1E-8;

   %% TolFun
   eps_f   = DefPar(options,'TolFun',optDef.TolFun);
   if eps_f ~= optDef.TolFun
      Prob.optParam.eps_f =  eps_f;
   end

   %% TolX
   eps_x   = DefPar(options,'TolX',optDef.TolX);
   if eps_x ~= optDef.TolX
      Prob.optParam.eps_x =  eps_x;
   end

   %% TolCon
   cTol   = DefPar(options,'TolCon',optDef.TolCon);
   if cTol ~= optDef.TolCon
      Prob.optParam.cTol =  cTol;
   end

   % MaxIter
   MaxIter = DefPar(options,'MaxIter',optDef.MaxIter);
   if MaxIter ~= optDef.MaxIter
      Prob.optParam.MaxIter = MaxIter;
   end
   % MaxFunEvals
   MaxFunc = DefPar(options,'MaxFunEvals',optDef.MaxFunEvals);
   if MaxFunc ~= optDef.MaxFunEvals
      Prob.optParam.MaxFunc = MaxFunc;
   end

   %% DiffMinChange
   DiffInt   = DefPar(options,'DiffMinChange',optDef.DiffMinChange);
   if DiffInt ~= optDef.DiffMinChange
      Prob.optParam.DiffInt =  DiffInt;
   end

   if isfield(options,'GradObj')
      if strcmpi('on',options.GradObj)
         Prob.NumDiff=0;
      elseif strcmpi('off',options.GradObj)
         Prob.NumDiff=1;
      end
   end
   if isfield(options,'GradConstr')
      if strcmpi('on',options.GradConstr)
         Prob.ConsDiff=0;
      elseif strcmpi('off',options.GradConstr)
         Prob.ConsDiff=1;
      end
   end
   if isfield(options,'Hessian')
      if strcmpi('on',options.Hessian)
         Prob.NumDiff=0;
      elseif strcmpi('off',options.Hessian)
         if Prob.NumDiff == 0
            Prob.NumDiff=-1;
         end
      end
   end

   Prob.HessPattern = DefPar(options,'HessPattern',Prob.HessPattern);
   Prob.HessMult = DefPar(options,'HessMult',[]);
   if ~isempty(Prob.HessMult)
      if Prob.Warning
       fprintf('\nWARNING! Inefficient to use option HessMult with Tomlab\n\n')
       fprintf('\nUse explicit sparse Hessian instead\n\n')
      end
   end

   if isfield(options,'LargeScale')
      if strcmpi('on',options.LargeScale)
         Prob.LargeScale=1;
      elseif strcmpi('off',options.LargeScale)
         Prob.LargeScale=0;
      end
   end
   if isfield(options,'LineSearchType')
      if strcmpi('cubicpoly',options.LineSearchType)
         Prob.LineParam.LineAlg=1;
      elseif strcmpi('quadcubic',options.LineSearchType)
         Prob.LineParam.LineAlg=0;
      end
   end
   if isfield(options,'Display')
      if strcmpi('iter',options.Display)
         Prob.optParam.IterPrint=1;
         PriLev=1;
      elseif strcmpi('final',options.Display)
         Prob.optParam.IterPrint=0;
         PriLev=1;
      elseif strcmpi('off',options.Display)
         Prob.optParam.IterPrint=0;
         PriLev=0;
      end
   end
   if isfield(options,'Diagnostics')
      if strcmpi('on',options.Diagnostics)
         Diagnostic=1;
      elseif strcmpi('off',options.Diagnostics)
         Diagnostic=0;
      end
   end
else
   Prob.HessMult = [];
end

if isa(Func,'cell')
   funArgIn  = 1;
   if length(Func) == 1
      Prob.OPTTB.f=Func{1};
      funArgOut = -1;
   elseif length(Func) == 2
      Prob.OPTTB.f=Func{1};
      Prob.OPTTB.g=Func{2};
      if isempty(Prob.OPTTB.g)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      else
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      end
   else
      Prob.OPTTB.f=Func{1};
      Prob.OPTTB.g=Func{2};
      Prob.OPTTB.H=Func{3};
      if isempty(Prob.OPTTB.g)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      elseif isempty(Prob.OPTTB.H)
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      else
         funArgOut = -3;
         Prob.NumDiff=0;
      end
   end
elseif isa(Func,'function_handle')
   Prob.OPTTB.f=Func;
   if M7
      Prob.OPTTB.M7 = 1;
      funArgIn  = abs(nargin(Func));
      funArgOut = abs(nargout(Func));
      if funArgOut == 1
         Prob.NumDiff = max(1,Prob.NumDiff);
      elseif funArgOut == 2
         Prob.NumDiff = min(-1,Prob.NumDiff);
      else
         Prob.NumDiff=0;
      end
   else
      SS = functions(Func);
      funArgIn  = abs(nargin(SS.function));
      funArgOut = abs(nargout(SS.function));
   end
%elseif isa(Func,'inline')
else
   Prob.OPTTB.f=Func;
   if strcmpi(Func,'nlp2_fgH')
      % Now Tomlab calls Tomlab through intermediate fmincon interface - weird
      % Must set proper NumDiff, in case derivatives not defined
      if isempty(Prob.FUNCSX.g)
         funArgOut = 1;
      elseif isempty(Prob.FUNCSX.H)
         funArgOut = 2;
      else
         funArgOut = 3;
      end
      funArgIn  = 3;
   else
      funArgOut = abs(nargout(Func));
      funArgIn  = abs(nargin(Func));
   end

   if funArgOut == 1
      Prob.NumDiff = max(1,Prob.NumDiff);
   elseif funArgOut == 2
      Prob.NumDiff = min(-1,Prob.NumDiff);
   else
      Prob.NumDiff=0;
   end
end

Prob.OPTTB.funArgOut=funArgOut;
Prob.OPTTB.funArgIn =funArgIn;

if isa(conFunc,'cell') & ~isempty(conFunc)
   conArgIn  = 1;
   if length(conFunc) == 1
      Prob.OPTTB.c=conFunc{1};
      conArgOut = -1;
   elseif length(conFunc) == 2
      Prob.OPTTB.c=conFunc{1};
      Prob.OPTTB.dc=conFunc{2};
      if isempty(Prob.OPTTB.dc)
         Prob.ConsDiff = max(1,Prob.ConsDiff);
         conArgOut = -1;
      else
         Prob.ConsDiff = min(-1,Prob.ConsDiff);
         conArgOut = -2;
      end
   else
      Prob.OPTTB.c=conFunc{1};
      Prob.OPTTB.dc=conFunc{2};
      Prob.OPTTB.d2c=conFunc{3};
      if isempty(Prob.OPTTB.dc)
         Prob.ConsDiff = max(1,Prob.ConsDiff);
         conArgOut = -1;
      elseif isempty(Prob.OPTTB.d2c)
         Prob.ConsDiff = min(-1,Prob.ConsDiff);
         conArgOut = -2;
      else
         conArgOut = -3;
         Prob.ConsDiff=0;
      end
   end
elseif isa(conFunc,'function_handle')
   Prob.OPTTB.c=conFunc;
   if M7
      Prob.OPTTB.M7Con = 1;
      conArgIn  = abs(nargin(conFunc));
      conArgOut = abs(nargout(conFunc));
      if conArgOut <= 2
         Prob.ConsDiff = max(1,Prob.ConsDiff);
      elseif conArgOut <= 4
         Prob.ConsDiff = min(-1,Prob.ConsDiff);
      else
         Prob.ConsDiff=0;
      end
   else
      SS = functions(conFunc);
      conArgIn  = abs(nargin(SS.function));
      conArgOut = abs(nargout(SS.function));
   end
%elseif isa(conFunc,'inline')
elseif ~isempty(conFunc)
   Prob.OPTTB.c=conFunc;
   if strcmpi(conFunc,'nlp2_cdceq')
      % Now Tomlab calls Tomlab through intermediate fmincon interface - weird
      % Must set proper ConsDiff, in case derivatives not defined
      if isempty(Prob.FUNCSX.dc)
         conArgOut = 1;
      elseif isempty(Prob.FUNCSX.d2c)
         conArgOut = 2;
      else
         conArgOut = 3;
      end
      conArgIn  = 3;
   else
      conArgOut = abs(nargout(conFunc));
      conArgIn  = abs(nargin(conFunc));
   end
   if conArgOut <= 2
      Prob.ConsDiff = max(1,Prob.ConsDiff);
   elseif conArgOut <= 4
      Prob.ConsDiff = min(-1,Prob.ConsDiff);
   else
      Prob.ConsDiff=0;
   end
else
   conArgOut = 0;
   conArgIn  = 0;
   Prob.ConsDiff=0;
end

Prob.OPTTB.conArgOut=conArgOut;
Prob.OPTTB.conArgIn =conArgIn;

% Must save r and J if the underlying problem is a least squares problem
r = Prob.FUNCS.r;
J = Prob.FUNCS.J;

if ~isempty(conFunc)
   if conArgOut < 0
      c=eval(conFunc);
      ceq = [];
   else
      if conArgIn > 1
         [c,ceq]=feval(conFunc,x_0,Prob.varargin{:});
      else
         [c,ceq]=feval(conFunc,x_0);
      end
   end
   c=c(:);
   ceq=ceq(:);
   mineq=length(c);
   meq=length(ceq);
   m=mineq+meq;
   Prob.c_L=[zeros(meq,1);-Inf*ones(mineq,1)];
   Prob.c_U=zeros(m,1);
   Prob = tomFiles(Prob,'optim_fgH','optim_g','optim_H','optim_cdc','optim_dc');
   Prob.mNonLin = m;
else
   Prob.ConsDiff = 0;
   Prob.mNonLin = 0;
   meq=0;mineq=0;
   Prob = tomFiles(Prob,'optim_fgH','optim_g','optim_H');
end

% Reset values
Prob.FUNCS.r=r;
Prob.FUNCS.J=J;
Prob = mkbound(Prob);
Prob.Solver.Alg  = DefPar(Prob.Solver,'Alg',0);

if strcmpi('nlpSolve',Solver)
   Alg=Prob.Solver.Alg;
   if isempty(Alg)
      Alg=0;
   end
   if Alg==0 & Prob.NumDiff~=0
      Prob.Solver.Alg=2;
   end
   if isfield(options,'HessUpdate')
      if strcmpi('bfgs',options.HessUpdate)
         Prob.Solver.Alg=2;
      elseif strcmpi('dfp',options.HessUpdate)
         Prob.Solver.Alg=2;
      end
   end

elseif strcmpi('conSolve',Solver)
   Alg=Prob.Solver.Alg;
   if isempty(Alg)
      Alg=0;
   end
   if Alg==0 & Prob.NumDiff~=0
      Prob.Solver.Alg=2;
   end
   if isfield(options,'HessUpdate')
      if strcmpi('bfgs',options.HessUpdate)
         Prob.Solver.Alg=2;
      elseif strcmpi('dfp',options.HessUpdate)
         Prob.Solver.Alg=2;
      end
   end
end

if Diagnostic
   Alg=Prob.Solver.Alg;
   if isempty(Alg)
      Alg=0;
   end
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf('   Diagnostic Information');
   fprintf('\n');
   if ~isempty(Prob.Name)
      fprintf('Problem %s',Prob.Name);
   end
   fprintf('\n');
   fprintf('Number of variables: %d\n',n);
   fprintf('Functions\n');
   fprintf(' Objective: TOMLAB optim_fgH is calling  %s',Func);
   fprintf('\n');
   fprintf(' Gradient:  TOMLAB optim_g ');
   if Prob.NumDiff==0
      fprintf(' using gradient in global variable NLP_g\n');
   else
      fprintf(' using finite-differencing method # %d\n',Prob.NumDiff);
   end
   if Prob.NumDiff==0
      fprintf(' Hessian:   TOMLAB optim_H is');
      fprintf(' using Hessian in global variable NLP_H\n');
   elseif Prob.NumDiff>0 & Alg==0
      fprintf(' Hessian:   TOMLAB optim_H is');
      fprintf(' using finite-differencing method # %d\n',Prob.NumDiff);
   elseif strcmpi(Solver,'conSolve') & Alg > 0
      fprintf(' Hessian:   TOMLAB is');
      fprintf(' using quasi-Newton update\n');
   elseif strcmpi(Solver,'nlpSolve') % Alg > 0
      fprintf(' Hessian:   TOMLAB is');
      fprintf(' using quasi-Newton update\n');
   end
   fprintf('\n');
   fprintf(' Number of lower bound constraints:              %d\n',...
           sum(~isinf(x_L)));
   fprintf(' Number of upper bound constraints:              %d\n',...
           sum(~isinf(x_U)));

   if isempty(Prob.b_U) | isempty(Prob.b_L)
      mEQ = 0;
   else
      mEQ=sum(~isinf(Prob.b_U) & Prob.b_L==Prob.b_U);
   end
   if isempty(Prob.b_U) 
      mU = 0;
   else
      mU =sum(~isinf(Prob.b_U))-mEQ;
   end
   if isempty(Prob.b_L) 
      mL = 0;
   else
      mL =sum(~isinf(Prob.b_L))-mEQ;
   end

   fprintf(' Number of linear equality constraints:          %d\n',mEQ);
   if mL > 0
      fprintf(' Number of linear lower bound constraints:       %d\n',mL);
   end
   fprintf(' Number of linear upper bound constraints:       %d\n',mU);

   if isempty(Prob.c_U) | isempty(Prob.c_L)
      mcEQ = 0;
   else
      mcEQ=sum(~isinf(Prob.c_U) & Prob.c_L==Prob.c_U);
   end
   if isempty(Prob.c_U) 
      mcU = 0;
   else
      mcU =sum(~isinf(Prob.c_U))-mcEQ;
   end
   if isempty(Prob.c_L) 
      mcL = 0;
   else
      mcL =sum(~isinf(Prob.c_L))-mcEQ;
   end

   fprintf(' Number of nonlinear equality constraints:       %d\n',mcEQ);
   if mL > 0
      fprintf(' Number of nonlinear lower bound constraints:    %d\n',mcL);
   end
   fprintf(' Number of nonlinear upper bound constraints:    %d\n',mcU);

   fprintf('\n');
   fprintf('Solver: %s',Solver);
   if isempty(Prob.Solver.Alg)
      fprintf('. Default algorithm.\n');
   else
      fprintf('. Alg # %d.\n',Prob.Solver.Alg);
   end
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf(' End diagnostic information\n');
   fprintf('\n');
end

Result   = tomRun(Solver,Prob,PriLev);

x        = x_0;          % Initialize x to x_0 to get correct size
x(:)     = Result.x_k;   % Column vector Result.x_k gets into x with size(x_0)
f_k      = Result.f_k;
g        = Result.g_k;
H        = Result.H_k;
ExitFlag = Result.ExitFlag;
Inform   = Result.Inform;

% Convert ExitFlag to OPTIM TB
if ExitFlag==0
   ExitFlag=1;
   if PriLev > 0
      fprintf('fmincon (%s',Solver);
      fprintf('): Optimization terminated successfully\n');
   end
else
   switch Inform
      case 101
        ExitFlag=0;
        if PriLev > 0
           fprintf('fmincon (%s',Solver);
           fprintf('): Too many iterations\n');
        end
      otherwise
        ExitFlag=-Inform;
        if PriLev > 0
           fprintf('fmincon (%s',Solver);
           fprintf('): Possible error.');
           ExitText = Result.ExitText;
           if ~isempty(ExitText)
              fprintf(' Solver information for Inform = %d:\n',Inform);
              for i = 1:size(ExitText,1)
                  fprintf('%s\n',ExitText(i,:));
              end
           else
              fprintf(' See solver help: Inform = %d\n',Inform);
           end
        end
   end
end

if nargout > 3
   Output.iterations  = Result.Iter;
   if isempty(Result.FuncEv)
      Output.funcCount   = n_f;
   else
      Output.funcCount   = Result.FuncEv;
   end
   Output.algorithm   = [ Solver ': ' Result.SolverAlgorithm];
   if isempty(alphaV)
      Output.stepsize   = 1;
   else
      Output.stepsize   = alphaV(length(alphaV));
   end
   Output.cgiterations  = [];
   Output.firstorderopt = [];
end

if nargout > 4
   [v_k, Zv, P] = LagMult(Result.Prob,Result);

   ix = find(P(1:n) > 0);
   Lambda.upper      = zeros(n,1);
   Lambda.upper(ix)  = v_k(ix);

   ix = find(P(1:n) == -1 | P(1:n) == 2);
   Lambda.lower      = zeros(n,1);
   Lambda.lower(ix)  = v_k(ix);

   Lambda.eqlin  = v_k(n+1:n+mAeq);
   
   n1=n+mAeq;

   Lambda.ineqlin      = zeros(mA,1);
   if mA > 0
      ix = find(P(n1+1:n1+mA) ~= 0);
      if ~isempty(ix)
         Lambda.ineqlin(ix)  = v_k(n1+ix);
      end
   end

   n1=n1+mA;
   Lambda.eqnonlin  = v_k(n1+1:n1+meq);
   n1=n1+meq;

   Lambda.ineqnonlin      = zeros(mineq,1);
   if mineq > 0
      ix = find(P(n1+1:n1+mineq) ~= 0);
      if ~isempty(ix)
         Lambda.ineqnonlin(ix)  = v_k(n1+ix);
      end
   end
end

% MODIFICATION LOG
%
% 011203 hkh  General solver handling with Prob.Solver.TOM
% 011204 hkh  Test is nonempty before setting value from options into Prob
% 011204 hkh  Save r and J for least squares problems
% 020416 hkh  Set Prob.x_0 before ProbCheck, then Prob.N is set.
% 020416 hkh  Revision regarding numerical derivatives
% 020416 hkh  Also test on options.LargeScale for solver selection
% 020929 hkh  Add TolCon, Prob.optParam.cTol
% 030114 hkh  Major revision to make fmincon work more efficient
% 030127 hkh  Handle case when x matrix
% 030127 hkh  Use global otxProb, and Prob.Solver.Tomlab for Tomlab solver name
% 030127 hkh  Use DefPar to set parameters
% 030128 hkh  Test for number of arguments when function_handle
% 030916 ango Fix missing [], line 499
% 031129 hkh  Empty x_0 gave wrong value of Set Prob.FUNCS.x, init x_0 first
% 031201 hkh  Change AutoDiff to fields ADObj and ADCons; error in options test
% 040124 hkh  Set Prob.A before call to ProbCheck. Set mNonLin field
% 040203 hkh  Use field FuncEv instead of global n_f
% 040203 hkh  Print ExitText if possible error (and PriLev > 0)
% 041222 hkh  Correct the defaults, change logic concerning defaults
% 041222 hkh  Call optParamDef with solver used, not with "quadprog"
% 041222 hkh  Handle new type of function handle in Matlab 7.x
% 050422 hkh  Test and set NumDiff/ConsDiff if calling Tomlab solver
% 060814 med  FUNCS used for callbacks instead
% 080414 hkh  Revise help, add text about choice of solver
% 080521 med  Defaults call added
% 080607 hkh  Use tomRun, not tomSolve
% 090813 med  Checks for x_0, x_L and x_U added
% 110821 hkh  Added comments for Lambda.eqnonlin och Lambda.ineqnonlin
% 110823 hkh  Also Lambda.lower in case of fixed variables


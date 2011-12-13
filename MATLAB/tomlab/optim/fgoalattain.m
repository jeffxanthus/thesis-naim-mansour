% FGOALATTAIN is the TOMLAB equivalent to fgoalattain in Optimization TB 
% The NLP solver actually used is selectable, see Prob.Solver.Tomlab below.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% TOMLAB fgoalattain solves the multi-objective goal attainment optimization
% problem on the OPTTB form:
%
%	min   {w: r(x) - w.*lam <= g)
%    x
%          x_L <=   x        <= x_U
%                   A   * x  <= b
%                   Aeq * x  == beq
%                   c(x)     <= 0
%                   ceq(x)   == 0
%
% The TOMLAB multi-objective goal attainment optimization is more general
% with respect to the constraints:
%
%	min   {w: r(x) - w.*lam <= g)
%    x
%          x_L <=   x   <= x_U
%          b_L <=   Ax  <= b_U
%          c_L <=  c(x) <= c_U
%
% where r is the m-dimensional function vector and
% x is a n-dimensional unknown parameter vector, with bounds x_L and x_U.
% A is a mA by n matrix with linear constraints, with bounds b_L and b_U
% c(x) is a m vector of nonlinear constraints, with bounds c_L and c_U
%
% The multi-objective goal attainment problem is solved as a minimax problem
% The OPTTB form is transformed to
%
%	min              gamma
%    x,gamma
%             x_L             <=   x   <= x_U
%             b_L = [-Inf ]   <=   A = [ A  ] * x       <= b_U = [ 0 ]
%                   [ beq ]   <=       [ Aeq]           <=       [beq]
%             c_L = [-Inf ]   <=   c = [ c(x)  ]        <= c_U = [ 0 ]
%                   [  0  ]   <=       [ ceq(x)]        <=       [ 0 ]
%                   [-Inf ]   <=       [ r(x)-w*gamma]  <=       [ g ]
%
% function [x, f_k, AttFact, ExitFlag, Output, Lambda] = fgoalattain(...
%    rFunc, x_0, GOAL, WEIGHT, A, b, Aeq, beq, x_L, x_U, conFunc, ...
%    options, Prob, varargin)
%
% INPUT: ( 4 arguments always needed )
% Func     Function that computes the m-dimensional function r(x), the set of
%          objectives.
%          Also optionally computes the Jacobian matrix corresponding to r(x)
% x_0      Initial values for the design variables
% GOAL     The set of design goals g. Dependent on the formulation, the
%          optimization will try to fulfill F <= GOAL, F == GOAL, or F >= GOAL
% WEIGHT   Set of weighting parameters w, which determine the relative under
%          or over achievement of the objectives.
%          1. Setting WEIGHT = abs(GOAL)  will try to make the objectives
%             less than the goals resulting in roughly the same
%             percentage under or over achievement of the goals.
%             Note: Use WEIGHT == 1 for GOAL values that are 0
%          2. Setting WEIGHT = -abs(GOAL) will try to make the objectives
%             greater then the goals resulting in roughly the same percentage
%             under- or over-achievement in the goals.
%             Note: use WEIGHT 1 for GOAL values that are 0
%          3. Setting WEIGHT(i)=0  indicates a hard constraint. i.e. F <= GOAL.
%             Such a hard constraint could equally well be added to the
%             nonlinear constraints
% A        Linear constraint matrix for inequalities
% b        Right hand side for linear inequality constraints
% Aeq      Linear constraint matrix for equalities
% beq      Right hand side for linear equality constraints
% x_L      Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% x_U      Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
% conFunc  Function that computes the constraints c(x) (and maybe the
%          constraint Jacobian dc(x))
%          If c(x) in R^m and x in R^n, then dc_k is an m-by-n matrix where
%          dc(i,j) is the partial derivative of c_i(x) with respect to x(j)
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX, TolFun, TolCon, DerivativeCheck,
%          GradObj, GradConstr,
%          (Hessian, HessPattern, HessUpdate, LineSearchType, )
%          MaxFunEvals, MaxIter, MeritFunction, GoalsExactAchieve,
%          (LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.)
%          Diagnostics, DiffMinChange and DiffMaxChange,
%
%          The field GoalsExactAchieve indicates the number of objectives for
%          which it is required for the objectives r(x) to equal the goals
%          GOAL).  Such objectives should be partitioned into the first few
%          elements of r.
%
%          Use the GradObj option to specify that Func may be called with two
%          two output arguments where the second, J, is the partial derivatives
%          of the function dr/dx, at the point x: [r,J] = feval(Func,x).
%          Use the GradConstr option to specify that conFunc may be called
%          with four output arguments:
%          [C,Ceq,GC,GCeq] = feval(conFunc,x) where GC is the partial
%          derivatives of the constraint vector of inequalities C and GCeq is
%          the partial derivatives of the constraint vector of equalities Ceq.
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
% Prob.Solver.Tomlab ='conSolve'; selects the Tomlab conSolve solver
% Prob.Solver.Tomlab ='snopt';    selects the Tomlab /SOL snopt solver
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
% Additional fields used in the Prob structure, beside standard fields:
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
% AttFact  Attainment factor at the solution x.
%          If AttFact is negative, the goals have been over-achieved;
%          if AttFact is positive, the goals have been under-achieved.
%
% ExitFlag exit condition of fgoalattain.
%      > 0 fgoalattain converged to a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 Errors, ExitFlag=-Inform, see the Inform parameter description
%
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.funcCount     Number of function evaluations
%   Output.algorithm     Type of algorithm used
%
%   Output.steplength    Length of last step
%   Output.cgiterations  Number of CG iterations (if used)
%   Output.firstorderopt The first-order optimality (if used)
%
% Lambda   Structure with Lagrange multipliers at the solution
%    Lambda.ineqlin   Lagrange multipliers for the linear inequalities A
%    Lambda.eqlin     Lagrange multipliers for the linear equalities Aeq
%    Lambda.lower     Lagrange multipliers for the lower bounds lb
%    Lambda.upper     Lagrange multipliers for the upper bounds ub
%
% g        The final gradient
% H        The final Hessian
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format)
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
%   Solver   TOMLAB solver used, e.g. conSolve or snopt
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Sept 29, 2002.     Last modified Jun 6, 2008.

function [x, f_k, AttFact, ExitFlag, Output, Lambda] = fgoalattain(...
   rFunc, x_0, GOAL, WEIGHT, A, b, Aeq, beq, x_L, x_U, conFunc, ...
   options, Prob, varargin)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

global otxProb alphaV n_f

if nargin < 13, Prob = [];
   if nargin < 12, options = [];
      if nargin < 11, conFunc = [];
         if nargin < 10, x_U = [];
            if nargin < 9, x_L = [];
               if nargin < 8, beq = [];
                  if nargin < 7, Aeq = [];
                     if nargin < 6, b = [];
                        if nargin < 5, A = []; 
            if nargin < 4 
  error('fgoalattain requires four input arguments Func, x_0, GOAL and WEIGHT');
end, end, end, end, end, end, end, end, end, end

% ------------Initialization----------------

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

z = version;
M7 = z(1) == '7';

solvType=checkType('cls');
Prob.probType = DefPar(Prob,'probType',solvType);

Prob.OPTTB.M7    = 0;
Prob.OPTTB.M7Con = 0;
Prob.OPTTB.x     = x_0;
Prob.x_0         = x_0(:);
Prob.x_L         = x_L(:);
Prob.x_U         = x_U(:);
n                = max([length(x_0(:)),length(x_L(:)),length(x_U(:))]);
Prob.N           = n;

Prob=ProbCheck(Prob,'fgoalattain',solvType);

mAeq=size(Aeq,1);
mA  =size(A,1);

b   = b(:);
beq = beq(:);

% Set default options

% TolX
Prob.optParam.eps_x=1E-6;

% TolFun
Prob.optParam.eps_f=1E-6;

% TolCon
Prob.optParam.cTol=1E-6;

% DerivativeCheck - off

% Diagnostics - off

Diagnostic=0;

% Hessian/gradient numerical approximation off
Prob.NumDiff=0; 
% Constraint Jacobian numerical approximation off
Prob.ConsDiff=0; 
% No Automatic Differentiation
Prob.ADObj=0;
Prob.ADCons=0;

% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

% LargeScale off as default
Prob.LargeScale=0;

% MaxFunEvals 100*numberOfVariables
% MaxIter 400
Prob.optParam.MaxIter=max(400,100*n);

% DiffMaxChange 1e-1 DiffMinChange 1e-8
Prob.optParam.DiffInt=1E-8;
%Prob.optParam.DiffGradMaxChange=1E-1;

% 'LineSearchType','quadcubic' default, but in TOMLAB we set Cubic as default
Prob.LineParam.LineAlg=1;

% LevenbergMarq    on

% PrecondBandWidth 0    TypicalX ones(numberOfVariables,1) 
% MaxPCGIter       max(1,floor(numberOfVariables/2))        TolPCG 0.1

if ~isempty(options)
   Prob.GoalsExact = DefPar(options,'GoalsExactAchieve',0);

   Prob.optParam.MaxIter = DefPar(options,'MaxIter',Prob.optParam.MaxIter);
   Prob.optParam.eps_f = DefPar(options,'TolFun',Prob.optParam.eps_f);
   Prob.optParam.eps_x = DefPar(options,'TolX',Prob.optParam.eps_x);
   Prob.optParam.cTol = DefPar(options,'TolCon',Prob.optParam.cTol);
   Prob.optParam.DiffInt = DefPar(options,'DiffGradMinChange',...
        Prob.optParam.DiffInt);

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
      fprintf('\nWARNING! Inefficient to use option HessMult with Tomlab\n\n')
      fprintf('\nUse explicit sparse Hessian instead\n\n')
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

Prob.LS.w = WEIGHT;
Prob.LS.g = GOAL;

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
Prob.mLin = mA+mAeq;

if isa(rFunc,'cell')
   funArgIn  = 1;
   if length(rFunc) == 1
      Prob.OPTTB.r=rFunc{1};
      funArgOut = -1;
   elseif length(rFunc) == 2
      Prob.OPTTB.r=rFunc{1};
      Prob.OPTTB.J=rFunc{2};
      if isempty(Prob.OPTTB.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      else
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      end
   else
      Prob.OPTTB.r=rFunc{1};
      Prob.OPTTB.J=rFunc{2};
      Prob.OPTTB.d2r=rFunc{3};
      if isempty(Prob.OPTTB.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      elseif isempty(Prob.OPTTB.d2r)
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      else
         funArgOut = -3;
         Prob.NumDiff=0;
      end
   end
elseif isa(rFunc,'function_handle')
   Prob.OPTTB.r=rFunc;
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
   Prob.OPTTB.r=rFunc;
   funArgOut = nargout(rFunc);
   funArgIn  = nargin(rFunc);
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
      conArgIn  = nargin(SS.function);
      conArgOut = nargout(SS.function);
   end
elseif ~isempty(conFunc)
   Prob.OPTTB.c=conFunc;
   conArgOut = nargout(conFunc);
   conArgIn  = nargin(conFunc);
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
   Prob = tomFiles(Prob,'optim_fgH','optim_g','optim_H', ...
         'optim_cdc','optim_dc', [], 'optim_rJ', 'optim_J');
   Prob.mNonLin = m;
   if m > 10 | n > 20
      Prob.ConsPattern = estConsPattern(Prob);
   else
      Prob.ConsPattern = ones(m,n);
   end
else
   Prob.ConsDiff = 0;
   Prob.mNonLin = 0;
   Prob = tomFiles(Prob,'optim_fgH','optim_g','optim_H', ...
   [], [], [], 'optim_rJ', 'optim_J');
end

Prob = mkbound(Prob);

Prob.Solver.Alg  = DefPar(Prob.Solver,'Alg',0);
Solver           = DefPar(Prob.Solver,'Tomlab',[]);

if isempty(Solver)
   Solver = GetSolver('con',n > 200 | Prob.LargeScale,0);
   Prob.Solver.Tomlab = Solver;
end

if strcmpi('conSolve',Solver)
   if Prob.NumDiff  == 6, Prob.NumDiff = 1; end
   if Prob.ConsDiff == 6, Prob.ConsDiff = 1; end
end

Prob.SolverInf = Solver;

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
   %fprintf(' Objective: TOMLAB optim_fgH is calling  %s',Func);
   fprintf(' Objective: TOMLAB optim_rJ,optim_J is calling  %s',rFunc);
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

Result   = goalSolve(Prob, 0, varargin{:});
PrintResult(Result,Prob.PriLev);

x        = x_0;
x(:)     = Result.x_k;
f_k      = Result.r_k;
ExitFlag = Result.ExitFlag;
Inform   = Result.Inform;
AttFact  = Result.f_k;

% Convert ExitFlag to OPTIM TB
if ExitFlag==0
   ExitFlag=1;
   if PriLev > 0
      fprintf('fgoalattain (%s',Solver);
      fprintf('): Optimization terminated successfully\n');
   end
else
   switch Inform
      case 101
        ExitFlag=0;
        if PriLev > 0
           fprintf('fgoalattain (%s',Solver);
           fprintf('): Too many iterations\n');
        end
      otherwise
        ExitFlag=-Inform;
        if PriLev > 0
           fprintf('fgoalattain (%s',Solver);
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

if nargout > 4
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

if nargout > 5
   v_k = Result.v_k;
   x_L = Result.Prob.x_L;
   x_U = Result.Prob.x_L;
   Lambda.upper      = zeros(n,1);
   ix                = find(x(:)==x_U);
   Lambda.upper(ix)  = v_k(ix);

   Lambda.lower      = zeros(n,1);
   ix                = find(x(:)==x_L);
   Lambda.lower(ix)  = v_k(ix);
   if Prob.mLin > 0
      b_L = Result.Prob.b_L;
      b_U = Result.Prob.b_U;
      ix                = find(b_L==b_U);
      Lambda.eqlin      = v_k(n+ix);
      ix                = find(b_L~=b_U);
      Lambda.ineqlin    = v_k(n+ix);
   else
      Lambda.eqlin      = [];
      Lambda.ineqlin    = [];
   end
   if Prob.mNonLin > 0
      c_L = Prob.c_L;
      c_U = Prob.c_U;
      ix                = find(c_L==c_U);
      Lambda.eqnonlin   = v_k(n+Prob.mLin+ix);
      ix                = find(c_L~=c_U);
      Lambda.ineqnonlin = v_k(n+Prob.mLin+ix);
   else
      Lambda.eqnonlin   = [];
      Lambda.ineqnonlin = [];
   end
end

% MODIFICATION LOG
%
% 020929 hkh Written, based on fmincon in Tomlab
% 030126 hkh Revised, based on numerous changes in fmincon
% 030128 hkh Further changes, similar to fmincon. Use DefPar, global struct.
% 030129 hkh Change name to optim_rJ from optim_r
% 031201 hkh Change AutoDiff to fields ADObj and ADCons; error in options test
% 040203 hkh Use field FuncEv instead of global n_f
% 040203 hkh Print ExitText if possible error (and PriLev > 0)
% 040215 hkh Call goalSolve with 0, not [], to avoid PrintResult
% 040526 hkh Set Prob.mLin and Prob.mNonLin
% 050422 hkh Handle new type of function handle in Matlab 7.x
% 050422 hkh Test and set NumDiff/ConsDiff if calling Tomlab solver
% 060817 hkh Use v_k from goalSolve instead of calling LagMult
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added
% quadprog is the TOMLAB equivalent to QUADPROG in Optimization TB
% The QP solver actually used is selectable, see Prob.SolverQP.
% If no active choice of solver is made, GetSolver will be called to select
% the best solver depending on the license used
%
% quadprog solves the quadratic programming problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x   <= b
%             x                                   Aeq*x == beq
%                                                 lb <= x <= ub
%
% function [x, fVal, ExitFlag, Output, Lambda, Result] = quadprog(H, f, ...
%           A, b, Aeq, beq, lb, ub, x0, options, Prob, varargin)
%
% INPUT: ( 3 arguments always needed )
%
% H        Symmetric matrix in the quadratic term in the objective function
% f        Linear term coefficients in the objective function
% A        Matrix of inequality constraint coefficients.
% b        Right hand side in inequality constraints
% Aeq      Matrix of equality constraint coefficients.
% beq      Right hand side in equality constraints
% lb       Lower bounds on the design values. -Inf == unbounded below.
%          Empty lb ==> -Inf on all variables
% ub       Upper bounds on the design values.  Inf == unbounded above.
%          Empty ub ==>  Inf on all variables
% x0       Set the starting point.
% options  Replaces the default optimization parameters
%          Fields used by Tomlab:
%                 Display, TolX, TolFun, Diagnostics, MaxIter, LargeScale
%          Fields not used by Tomlab:
%                 MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%          HessMult: If nonempty, Tomlab makes i=1:dim(x) calls
%             feval(HessMult,H,e_i,varargin),   to obtain the H matrix.
%          e_i= (....,1, ....), 1 in the ith position, 0:s otherwise
%          The Tomlab solvers cannot utilize the H*x computation in HessMult
%
% Prob     The TOMLAB problem input structure, or the first extra user
%          input argument
%
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
% Additional input as fields in Prob:
%
% Prob.SolverQP   Name of Tomlab QP solver to use
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.SolverQP = 'qpopt';
%
% Note! If the QP problem is convex, suitable TOMLAB solvers are
%       cplex (large-scale), sqopt (large-scale), qp-minos (large-scale),
%       qpopt(medium size) or qpSolve(small size).
%       If the QP problem is nonconvex, suitable TOMLAB solvers are
%       snopt (large-scale), knitro (large-scale), minos (large-scale),
%       qpopt(medium size) or qpSolve(small size).
%
% OUTPUT:
%
% x        Optimal design parameters
% fVal     Optimal design parameters
% ExitFlag exit condition of quadprog.
%      > 0 quadprog converged with a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 No convergence
%      - 2 Unbounded feasible region.
%      - 3 Rank problems
%      - 4 Illegal x0 found
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.algorithm     Type of algorithm used
%   Output.cgiterations  Number of CG iterations (LargeScale on)
%   Output.firstorderopt First order optimality  (LargeScale on)
% Lambda   Structure with Lagrange multipliers at the solution
%   Lambda.ineqlin   Lagrange multipliers for the linear inequalities A
%   Lambda.eqlin     Lagrange multipliers for the linear equalities Aeq
%   Lambda.lower     Lagrange multipliers for the lower bounds lb
%   Lambda.upper     Lagrange multipliers for the upper bounds ub
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format). Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%            == 1  => Maximal number of iterations reached. No bfs found.
%            == 2  => Unbounded feasible region.
%            == 3  => Rank problems
%            == 4  => No feasible point found with Phase 1 algorithm
%   Inform   If ExitFlag > 0, Inform=ExitFlag, otherwise Inform show type
%            of convergence:
%            0 = Unconstrained solution
%            1 = lambda >= 0.
%            2 = lambda >= 0. No 2nd order Lagrange mult. estimate available
%            3 = lambda and 2nd order Lagrange mult. positive, problem is
%                not negative definite.
%            4 = Negative definite problem. 2nd order Lagrange mult. positive,
%                but releasing variables leads to same working set.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   f_k      Function value 0.5*x'*F*x+c'*x
%   g_k      Gradient F*x+c
%   H_k      Hessian F (constant)
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   Solver   The TOMLAB solver used, e.g. cplex, snopt, sqopt or qpSolve
%   SolverAlgorithm  Description of method used

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written July 3, 1999.   Last modified Aug 13, 2009.

function [x, fVal, ExitFlag, Output, Lambda, Result] = quadprog(H, f, ...
          A, b, Aeq, beq, lb, ub, x0, options, Prob, varargin)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

global otxProb

if nargin < 11, Prob = [];
   if nargin < 10, options = [];
      if nargin < 9, x0 = []; 
         if nargin < 8, ub = []; 
            if nargin < 7, lb = []; 
               if nargin < 6, beq = [];
                  if nargin < 5, Aeq = [];
                     if nargin < 4, b = [];
                        if nargin < 3, A = [];
                           if nargin < 2
                              error('quadprog requires two parameters H,f');
end, end, end, end, end, end, end, end, end, end

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

HessMult = DefPar(options,'HessMult',[]);

if ~isempty(HessMult)
   n = max([length(lb),length(f),length(ub)]);
   y = zeros(n,1);
   Z = H;
   H = sparse(n,n);
   % Sick way to obtain the H matrix - works, but slow 
   for i=1:n
       y(i) = 1;
       H(:,i) = feval(HessMult, Z, y, Prob.varargin{:});
       y(i) = 0;
   end
   H = 0.5*(H+H');
else
   n = max([size(H),length(f)]);
   % Check the symmetry
   if ~isempty(H)
      if any(any(H-H' ~= 0))
         if Prob.Warning
            warning('The Hessian is not symmetric! Using H=(H+H'')/2 instead');
         end
         H = 0.5*(H+H');
      end
   end
end

Prob.QP.F= H;
Prob.QP.c= f(:);    % cost vector

solvType=checkType('qp');
if isempty(H)
   Prob.probType = checkType('lp');
   Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');
else
   Prob.probType = solvType;
   Prob = tomFiles(Prob, 'qp_f', 'qp_g', 'qp_H');
end

n1  = size(A,1);
n2  = size(Aeq,1);
if isempty(Aeq)
   Prob.A   = A;
elseif isempty(A)
   Prob.A   = Aeq;
else
   Prob.A   = [A;Aeq];
end

Prob.b_U = cat(1,b(:),beq(:));

if isempty(beq)
   Prob.b_L = -Inf*ones(n1,1);
else
   Prob.b_L = [-Inf*ones(n1,1);beq(:)];
end

if length(x0) < n
    Prob.x_0 = zeros(n,1);
else
    Prob.x_0 = x0(:);
end

if length(lb) < n
    Prob.x_L = -inf*ones(n,1);
else
    Prob.x_L = lb(:);
end

if length(ub) < n
    Prob.x_U = inf*ones(n,1);
else
    Prob.x_U = ub(:);
end

Prob.N = n;
Prob.mLin    = n1+n2;
Prob.mNonLin = 0;

if isempty(Prob.SolverQP)
   Prob.SolverQP = GetSolver('qp',n > 200,1);
   %if strcmpi(Prob.SolverQP,'qpopt')
   %end
end

Solver=Prob.SolverQP;

% Set default options
% LargeScale on by default
Prob.LargeScale=1;

% Display final is default
Prob.optParam.IterPrint=0;

%Prob=ProbCheck(Prob,'quadprog',solvType);
Prob=ProbCheck(Prob,Solver,solvType);

Diagnostic=0;

% Opt tbx defaults
% Display          final
% Diagnostics      off
% LargeScale       on
% Opt tbx defaults (not used)
% HessMult         []
% MaxIter          200
% PrecondBandWidth 0    
% TypicalX         ones(numberOfVariables,1) 
% TolPCG           0.1
% MaxPCGIter       max(1,floor(numberOfVariables/2))

if ~isempty(options)
   % Opt tbx defaults - test if changed
   optDef.TolFun  = 100*eps;
   optDef.TolX    = 100*eps;
   optDef.MaxIter = 200;

   %% TolFun
   %Prob.optParam.eps_f=100*eps;
   eps_f   = DefPar(options,'TolFun',optDef.TolFun);
   if eps_f ~= optDef.TolFun
      Prob.optParam.eps_f =  eps_f;
   end

   %% TolX
   %Prob.optParam.eps_x=100*eps;
   eps_x   = DefPar(options,'TolX',optDef.TolX);
   if eps_x ~= optDef.TolX
      Prob.optParam.eps_x =  eps_x;
   end

   % MaxIter 200 in opt tbx
   % Prob.optParam.MaxIter=;
   MaxIter = DefPar(options,'MaxIter',optDef.MaxIter);
   if MaxIter ~= optDef.MaxIter
      % Only allow the user to set MaxIter >= 2000
      Prob.optParam.MaxIter = max(2000,MaxIter);
   end

   if isfield(options,'LargeScale')
      if strcmpi('on',options.LargeScale)
         Prob.LargeScale=1;
      else
         Prob.LargeScale=0;
      end
   end
   if isfield(options,'Display')
      if strcmpi('iter',options.Display)
         Prob.optParam.IterPrint=1;
      elseif strcmpi('final',options.Display)
         Prob.optParam.IterPrint=0;
      elseif strcmpi('off',options.Display)
         Prob.optParam.IterPrint=0;
      end
   end
   if isfield(options,'Diagnostics')
      if strcmpi('on',options.Diagnostics)
         Diagnostic=1;
      elseif strcmpi('off',options.Diagnostics)
         Diagnostic=0;
      end
   end
end

Prob = mkbound(Prob);

if Diagnostic
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf('   Diagnostic Information');
   fprintf('\n');
   if ~isempty(Prob.Name)
      fprintf('Problem %s',Prob.Name);
   end
   fprintf('\n');
   fprintf('Number of variables: %d\n',n);
   fprintf('\n');
   fprintf(' Number of lower bound constraints:              %d\n',...
           sum(~isinf(Prob.x_L)));
   fprintf(' Number of upper bound constraints:              %d\n',...
           sum(~isinf(Prob.x_U)));
   mEQ=sum(~isinf(Prob.b_U) & Prob.b_L==Prob.b_U);
   mU =sum(~isinf(Prob.b_U))-mEQ;
   mL =sum(~isinf(Prob.b_L))-mEQ;
   fprintf(' Number of linear equality constraints:          %d\n',mEQ);
   if mL > 0
      fprintf(' Number of linear lower bound constraints:       %d\n',mL);
   end
   fprintf(' Number of linear upper bound constraints:       %d\n',mU);
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

if strcmpi('qpSolve',Solver)
   Result   = qpSolve(Prob);
else
   Result   = tomRun(Solver,Prob,0);
end

PrintResult(Result,Prob.PriLev)

x        = Result.x_k;
fVal     = Result.f_k;
ExitFlag = Result.ExitFlag;

% Convert ExitFlag to OPTIM TB 
switch ExitFlag
case 0 
     ExitFlag=1;
case 1 
     ExitFlag=0;
case 2 
     ExitFlag=-2;
     fprintf('quadprog (%s): Unbounded feasible region\n',Solver);
case 3 
     ExitFlag=-3;
     fprintf('quadprog (%s): Rank problems\n',Solver);
case 4 
     ExitFlag=-4;
     fprintf('quadprog (%s): The problem is infeasible\n',Solver);
end

Output.iterations  = Result.Iter;
Output.algorithm   = [Solver ': ' Result.SolverAlgorithm];
if Prob.LargeScale
   Output.cgiterations = [];
   Output.firstorderopt= [];
end

if nargout > 4
   [v_k, Zv, P] = LagMult(Result.Prob,Result);

   if n1 > 0
   Lambda.ineqlin = v_k(n+1:n+n1);
   else
      Lambda.ineqlin = [];
   end
   if n2 > 0
      Lambda.eqlin   = -v_k(n+n1+1:n+n1+n2);
   else
      Lambda.eqlin   = [];
   end

   ix = find(P(1:n) > 0);
   Lambda.upper      = zeros(n,1);
   Lambda.upper(ix)  = v_k(ix);

   ix = find(P(1:n) == -1);
   Lambda.lower      = zeros(n,1);
   Lambda.lower(ix)  = v_k(ix);
end

% MODIFICATION LOG
%
% 011203 hkh Revision. Use Prob.SolverQP, makes the QP solver selectable
% 020821 hkh Must test if options.eps_x and others are empty
% 030128 hkh Add Prob given as global structure otxProb
% 030128 hkh Generate H if HessMult option is set
% 040303 hkh Define m-files for LP or QP, if general NLP solvers are used
% 041221 hkh Correct the defaults, change logic concerning defaults
% 041221 hkh Call optParamDef with solver used, not with "quadprog"
% 041221 hkh norm(abs(H-H'),inf) > eps scaling dependent, use any(any(F-F'))~=0
% 080414 hkh Assume default QP Problem is convex when calling GetSolver
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added
% 090127 med Return code infeasible problems corrected
% 090813 med x_0, x_L and x_U checks added
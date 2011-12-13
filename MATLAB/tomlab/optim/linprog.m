% linprog is the TOMLAB equivalent to LINPROG in Optimization TB
% The LP solver actually used is selectable, see Prob.SolverLP.
% If no active choice of solver is made, GetSolver will be called to select
% the best solver depending on the license used
%
% linprog solves the linear programming problem:
%
%            min f'*x    subject to:   A*x   <= b
%             x                        Aeq*x == beq
%                                      lb <= x <= ub
%
% function [x, fVal, ExitFlag, Output, Lambda, Result] = linprog(f, A, b, ...
%           Aeq, beq, lb, ub, x0, options, Prob)
%
% INPUT: ( 3 arguments always needed )
%
% f        Objective function coefficients
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
%          Fields used: Display, TolFun, Diagnostics, MaxIter,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%
% Prob     The TOMLAB problem input structure (default empty)
%
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
% Additional input as fields in Prob:
%
% Prob.SolverLP   Name of LP solver to use
%
% Suitable large-scale TOMLAB solvers are cplex, lp-minos, sqopt, conopt
%
% Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.SolverLP = 'sqopt';
%
% NOTE! linprog only checks for the global otxProb is input Prob is empty
%
% OUTPUT:
%
% x        Optimal design parameters
% fVal     Optimal design parameters
% ExitFlag exit condition of linprog.
%      > 0 linprog converged with a solution X.
%        0 Reached the maximum number of iterations without convergence
%      - 2 Unbounded feasible region.
%      - 3 Rank problems
%      - 4 Illegal x0 found
%      - 5 No feasible point found
% Output   Structure. Fields:
%   Output.iterations   Number of iterations
%   Output.algorithm    Type of algorithm used
%   Output.cgiterations Number of CG iterations (LargeScale on)
% Lambda   Structure with Lagrange multipliers at the solution
%    Lambda.ineqlin   Lagrange multipliers for the linear inequalities A
%    Lambda.eqlin     Lagrange multipliers for the linear equalities Aeq
%    Lambda.lower     Lagrange multipliers for the lower bounds lb
%    Lambda.upper     Lagrange multipliers for the upper bounds ub
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format). Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   f_k      Function value 0.5*x'*F*x+c'*x
%   g_k      Gradient c
%   Solver   Tomlab solver, e.g. cplex, sqopt or lp-minos
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written July 2, 1999.   Last modified Aug 13, 2009.

function [x, fVal, ExitFlag, Output, Lambda, Result] = linprog(f, A, b, ...
          Aeq, beq, lb, ub, x0, options, Prob)

global otxProb

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end
      
if nargin < 10, Prob = [];
   if nargin < 9, options = [];
      if nargin < 8, x0 = []; 
         if nargin < 7, ub = []; 
            if nargin < 6, lb = []; 
               if nargin < 5, beq = [];
                  if nargin < 4, Aeq = [];
                     if nargin < 3
                        error('linprog requires three parameters f,A,b');
end, end, end, end, end, end, end, end

if isempty(Prob)
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
end
Prob.probType=checkType('lp');

Prob.QP.c= f(:);    % cost vector
n        = length(f(:));
Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');
m1  = size(A,1);
m2  = size(Aeq,1);
if isempty(Aeq)
   Prob.A   = A;
elseif isempty(A)
   Prob.A   = Aeq;
else
   Prob.A   = [A;Aeq];
end
Prob.mLin = size(Prob.A,1);
Prob.mNonLin = 0;
Prob.b_U = cat(1,b(:),beq(:));

if isempty(beq)
   Prob.b_L = -Inf*ones(m1,1);
else
   Prob.b_L = [-Inf*ones(m1,1);beq(:)];
end

Prob.QP.F=[];

% LargeScale on
Prob.LargeScale=1;
% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;
N = max([length(lb),length(ub),length(f)]);
Prob.N = N;
if length(lb) < N
    Prob.x_L = -inf*ones(N,1);
else
    Prob.x_L = lb(:);
end
if length(ub) < N
    Prob.x_U = inf*ones(N,1);
else
    Prob.x_U = ub(:);
end
if length(x0) < N
    Prob.x_0 = zeros(N,1);
else
    Prob.x_0 = x0(:);
end

if isempty(Prob.SolverLP)
   Prob.SolverLP = GetSolver('lp',n > 200,0);
end

Solver=Prob.SolverLP;
Prob = ProbCheck(Prob,Solver,8);

% Set default options as in opt tbx
% Opt tbx defaults
% Display          final
% Diagnostics      off
% LargeScale       on
% MaxIter          [], (85 in v2.x)
%                  Actual default is either 100*Number of Variables or 85
% TolFun           [], (1E-8 in v2.x)
%                  Actual default is either 1E-6 or 1E-8
% Opt tbx defaults (not used)
% Simplex          off

% Diagnostics - off
Diagnostic=0;

if ~isempty(options)
   z = version;
   if z(1) == '6'
      % Matlab 6.x, opt tbx 2.x defaults
      % Opt tbx defaults - test if changed
      optDef.TolFun  = 1E-8;
      optDef.MaxIter = 85;
   else
      % Matlab 7.x, opt tbx 3.x defaults
      optDef.MaxIter = [];
      optDef.TolFun  = [];
   end

   %% TolFun
   %Prob.optParam.eps_f=1E-8;
   eps_f   = DefPar(options,'TolFun',optDef.TolFun);
   if eps_f ~= optDef.TolFun
      Prob.optParam.eps_f =  eps_f;
   end

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
elseif strcmpi('lpSimplex',Solver)
   Result   = lpSimplex(Prob);
elseif strcmpi('minos',Solver)
   Prob.SOL.optPar(30)=100000;
   if Prob.SOL.optPar(48) <= 50
      Prob.SOL.optPar(48)=500;
   end
   Result   = tomRun(Solver,Prob);
elseif strcmpi('lp-minos',Solver)
   Prob.SOL.optPar(30)=100000;
   Result   = tomRun(Solver,Prob);
elseif strcmpi('qp-minos',Solver)
   Prob.SOL.optPar(30)=100000;
   Result   = tomRun(Solver,Prob);
else
   Result   = tomRun(Solver,Prob);
end

% Result is last output variable.
% It is then easy to call PrintResult with this Result output
%
% Prob.PriLev=2
PrintResult(Result,Prob.PriLev);

x        = Result.x_k;
fVal     = Result.f_k;
ExitFlag = Result.ExitFlag;

% Convert ExitFlag to OPTIM TB 
switch ExitFlag
case 0 
     ExitFlag=1;
     if PriLev > 0
        fprintf('linprog (%s',Solver);
        fprintf('): Optimization terminated successfully\n');
     end
case 1 
     ExitFlag=0;
case 2 
     ExitFlag=-2;
     if PriLev > 0
        fprintf('linprog (%s',Solver);
        fprintf('): Unbounded feasible region\n');
     end
case 3 
     ExitFlag=-3;
     if PriLev > 0
        fprintf('linprog (%s',Solver);
        fprintf('): Rank problems\n');
     end
case 4 
     ExitFlag=-4;
     if PriLev > 0
        fprintf('linprog (%s',Solver);
        fprintf('): The problem is infeasible\n');
     end
case 5 
     ExitFlag=-5;
     if PriLev > 0
        fprintf('linprog (%s',Solver);
        fprintf('): No feasible point found with Phase 1 simplex\n');
     end
end

if strcmpi('minos',Solver) | strcmpi('lp-minos',Solver) ...
                           | strcmpi('qp-minos',Solver)
   Output.iterations  = Result.MinorIter;
else
   Output.iterations  = Result.Iter;
end
Output.algorithm   = [ Solver ': ' Result.SolverAlgorithm];
if Prob.LargeScale
   Output.cgiterations= 0;
end

if nargout > 4
   v_k = Result.v_k;

   if m1 > 0
      Lambda.ineqlin = v_k(n+1:n+m1);
   else
      Lambda.ineqlin = [];
   end
   if m2 > 0
      Lambda.eqlin   = -v_k(n+m1+1:n+m1+m2);
   else
      Lambda.eqlin   = [];
   end

   ix = find(Result.Prob.x_U(1:n) == x(1:n));
   Lambda.upper      = zeros(n,1);
   Lambda.upper(ix)  = v_k(ix);

   ix = find(Result.Prob.x_L(1:n) == x(1:n));
   Lambda.lower      = zeros(n,1);
   Lambda.lower(ix)  = v_k(ix);
end

% MODIFICATION LOG
%
% 990729 hkh Add usage of qpSolve
% 990812 hkh Use Lagrange multipliers from the output of lpSimplex.
% 011203 hkh Revision. Use Prob.SolverLP, makes the LP solver selectable
% 030128 hkh Add option Prob given as global structure otxProb
% 030502 hkh Increase Superbasics limit as default for MINOS
% 040303 hkh Set fields in Prob.FUNCS (f,g,H)
% 041221 hkh Correct the defaults, change logic concerning defaults
% 041221 hkh Call optParamDef with solver used, not with "linprog"
% 041221 hkh For minos,lp-minos,qp-minos, report MinorIter, not Iter
% 060819 med Removed Output print
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added
% 080917 med Return code for ExitFlag 4 corrected
% 090813 med x_0, x_L and x_U checks added
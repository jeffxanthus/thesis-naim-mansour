% bintprog is the TOMLAB equivalent to BINTPROG in Optimization TB 
% The MIP solver actually used is selectable, see Prob.SolverMIP.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% bintprog solves binary linear programming problems:
%
%            min f'*x    subject to:   A*x   <= b
%             x                        Aeq*x == beq
%
%            where the elements of x are binary
%            integers, i.e., 0's or 1's.
%
% function [x, fVal, ExitFlag, Output, Result] = bintprog(f, A, b, ...
%           Aeq, beq, x0, options, Prob)
%
% INPUT: ( 1 argument always needed )
%
% f        Objective function coefficients
% A        Matrix of inequality constraint coefficients
% b        Right hand side in inequality constraints
% Aeq      Matrix of equality constraint coefficients
% beq      Right hand side in equality constraints
% x0       Set the starting point.
% options  Replaces the default optimization parameters
%          Fields used: Display, TolFun, Diagnostics, MaxIter.
%
% Prob     The TOMLAB problem input structure (default empty)
%
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
%          Prob.PriLev is used if printing from TOMLAB is needed.
%
% Additional input as fields in Prob:
%
% Prob.SolverMIP   Name of MILP(MIP) solver to use
%
% Another way to input the TOMLAB problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty TOMLAB Prob structure
%          "set fields in otxProb", e.g. otxProb.SolverMIP = 'cplex';
%
% NOTE! bintprog only checks for the global otxProb is input Prob is empty
%
% OUTPUT:
%
% x        Optimal design parameters
% fVal     Optimal value
% ExitFlag exit condition of bintprog.
%        1 bintprog converged to a solution X.
%        0 Maximum number of iterations exceeded.
%      - 2 Problem is infeasible.
%      - 4 MaxNodes reached without converging.
%      - 5 MaxTime reached without converging.
%      - 6 Number of iterations performed at a node to solve the LP-relaxation
%          problem exceeded MaxRLPIter reached without converging.
%
% Output   Structure. Fields:
%   Output.iterations   Number of iterations
%   Output.algorithm    Type of algorithm used
%
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format). Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   QP.B     B  Optimal set. B(i)==0, include variable x(i) in basic set.
%            sum(B==1)==length(b)  holds.
%   f_k      Function value
%   g_k      Gradient c
%   Solver   Any TOMLAB solver, e.g. CPLEX
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2010 by Tomlab Optimization Inc., $Release: 7.5.0$
% Written Aug 29, 2005.   Last modified Sep 23, 2010.

function [x, fVal, ExitFlag, Output, Result] = bintprog(f, A, b, ...
          Aeq, beq, x0, options, Prob)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

global otxProb       
       
if nargin < 8, Prob = [];
   if nargin < 7, options = [];
      if nargin < 6, x0 = []; 
         if nargin < 5, beq = []; 
            if nargin < 4, Aeq = []; 
               if nargin < 3, b = [];
                  if nargin < 2, A = [];
                     if nargin < 1
                        error('bintprog requires one parameter f');
end, end, end, end, end, end, end, end

if isempty(Prob)
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
end

Prob.probType=checkType('mip');
Prob.QP.B = [];
Prob.QP.c = f(:);    % cost vector
n         = length(f(:));
% Define n variables as integer valued
Prob.MIP = struct('IntVars',double(1:n),...
         'VarWeight',[],'KNAPSACK',[],'fIP',[], 'xIP',[],'PI',[],...
         'SC',[],'SI',[],'sos1',[],'sos2',[],'xpcontrol',[],'callback',[]);
Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');
m1 = size(A,1);
if isempty(Aeq)
   Prob.A   = A;
elseif isempty(A)
   Prob.A   = Aeq;
else
   Prob.A   = [A;Aeq];
end
Prob.b_U = [b(:);beq(:)];
Prob.mLin = size(Prob.A,1);
Prob.mNonLin = 0;
if isempty(beq)
   Prob.b_L = -Inf*ones(m1,1);
else
   Prob.b_L = [-Inf*ones(m1,1);beq(:)];
end
Prob.x_0 = x0(:);
Prob.x_L = zeros(n,1);
Prob.x_U = ones(n,1);
Prob.N   = n;
Prob.QP.F=[];
Prob.HessPattern=sparse(n,n);

% LargeScale on
Prob.LargeScale=1;
% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

if ~isfield(Prob, 'SolverMIP')
   Prob.SolverMIP = [];
end
if isempty(Prob.SolverMIP)
   Prob.SolverMIP = GetSolver('mip');
end
Solver=Prob.SolverMIP;
Prob = ProbCheck(Prob,Solver,7);

% Set default options as in opt tbx
% Opt tbx defaults
% Display: 'final'                       DONE
% MaxIter: '100000*numberofvariables'    DONE
% TolFun: 1.0000e-003                    DONE
% BranchStrategy: 'maxinfeas'
% Diagnostics: 'off'                     DONE
% MaxNodes: '1000*numberofvariables'
% MaxRLPIter: '100*numberofvariables'
% MaxTime: 7200
% NodeSearchStrategy: 'bn'
% TolRLPFun: 1.0000e-006
% TolXInteger: 1.0000e-008

Diagnostic=0;

if ~isempty(options)
   optDef.TolFun  = 1E-3;
   optDef.MaxIter = 100000*n;

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
      fprintf(' Problem: %s',Prob.Name);
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
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf(' End diagnostic information\n');
   fprintf('\n');
end

Result   = tomRun(Solver,Prob);
x        = Result.x_k;
if nargin > 1
   fVal     = Result.f_k;
end
ExitFlag = Result.ExitFlag;

% Convert ExitFlag to OPTIM TB
switch ExitFlag
   case 0
      ExitFlag=1;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): bintprog converged to a solution X.\n');
      end
   case 1
      ExitFlag=0;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): Maximum number of iterations exceeded.\n');
      end
   case 2
      ExitFlag=-2;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): Problem is infeasible.\n');
      end
   case 3
      ExitFlag=-3;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): Rank problems\n');
      end
   case 4
      ExitFlag=-2;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): Problem is infeasible.\n');
      end
   case 5
      ExitFlag=-6;
      if PriLev > 0
         fprintf('bintprog (%s',Solver);
         fprintf('): Number of iterations performed at a node to solve the LP-relaxation\n');
         fprintf('problem exceeded MaxRLPIter reached without converging.');
      end
end

if nargin > 3
   Output.iterations  = Result.Iter;
   Output.algorithm   = [ Solver ': ' Result.SolverAlgorithm];
end

% MODIFICATION LOG
%
% 050829 med Created
% 050830 med Added isfield check for SolverMIP
% 050830 med General mlint fixes
% 070221 hkh Change to new IntVars format
% 080414 hkh Revise comments, add text about choice of solver
% 080521 med Defaults call added
% 080909 med Updated Prob.N to be correct
% 100923 hkh Define Prob.QP.B, because check is deleted from ProbCheck
% 100923 hkh Define a full proper Prob.MIP

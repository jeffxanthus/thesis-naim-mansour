% lsqlin is the TOMLAB equivalent to LSQLIN in Optimization TB
% The LLS/QP solver actually used is selectable, see Prob.SolverQP.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% TOMLAB lsqlin solves linear least squares problems on the form:
%
%	min       0.5 * sum ||C*x -d||^2 = f(x)
%        x
%                 x_L <=        x  <= x_U
%                         A   * x  <= b
%                         Aeq * x  <= beq
%
% x is a n-dimensional unknown parameter vector, with bounds x_L and x_U.
%
% function [x, f_k, r_k, ExitFlag, Output, Lambda] = lsqlin(...
%           C, d, A, b, Aeq, beq, x_L, x_U, x_0, options, Prob, varargin)
%
% INPUT: ( 2 arguments always needed )
%
% C        Matrix C in C*x-d.
% d        Data vector in C*x-d
% A        Linear constraint matrix for inequalities
% b        Right hand side for linear inequality constraints
% Aeq      Linear constraint matrix for equalities
% beq      Right hand side for linear equality constraints
% x_L      Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% x_U      Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
% x_0      Starting value for the design variables
% options  Replaces the default optimization parameters
%          Fields used: Display, TolFun, Diagnostics, MaxIter,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%          JacobMult: If nonempty, Tomlab makes i=1:dim(x) calls
%             feval(JacobMult,J,e_i,-1,varargin),   to obtain the J matrix.
%          e_i= (....,1, ....), 1 in the ith position, 0:s otherwise
%          The Tomlab solvers cannot utilize the J*x computation in JacobMult
%          (except for Tlsqr)
%
% Prob     The TOMLAB problem input structure, or the first extra user
%          input argument
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
%          If Prob.SolverQP  If set, this solver is used instead of the default solver
%
%          Weighting, using input weightType and weightY (and data vector y) is
%          normally computed in llsAssign.
%          Running lsqlin, the user must apply the same weighting directly to C and d
%
%          If set, the following fields are used in Prob:
%
%              Prob.SOL (used by SOL solvers)
%              Prob.optParam
%              Prob.LargeScale
%              Prob.LineParam.LineAlg (only used by clsSolve)
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. 
%                otxProb.SolverQP = 'lssol';
%          makes lsqlin call LSSOL to solve the linear least squares problem
%
% OUTPUT:
%
% x        Optimal design parameters
% f_k      Optimal residual sum of squares, sum {r(x).^2} (Note! no 0.5)
% r_k      The optimal residual vector
% ExitFlag exit condition of lsqlin.
%      > 0 lsqlin converged to a solution X.
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
%   Lambda.lower     Lagrange multipliers for the lower bounds lb
%   Lambda.upper     Lagrange multipliers for the upper bounds ub

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.7.0$
% Written July 29, 1999.  Last modified July 22, 2011.

function [x, f_k, r_k, ExitFlag, Output, Lambda] = lsqlin(...
          C, d, A, b, Aeq, beq, x_L, x_U, x_0, options, Prob,varargin)

global otxProb      

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end
   
if nargin < 11, Prob = [];
   if nargin < 10, options = [];
      if nargin < 9, x_0 = []; 
         if nargin < 8, x_U = []; 
            if nargin < 7, x_L = []; 
               if nargin < 6, beq = []; 
                  if nargin < 5, Aeq = []; 
                     if nargin < 4, b = []; 
                        if nargin < 3, A = []; 
            if nargin < 2 
	       error('lsqlin requires two input arguments C and d');
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

JacobMult = DefPar(options,'JacobMult',[]);

if ~isempty(JacobMult)
   n = max([length(x_L),length(x_L),length(x_0),size(A,2)]);
   y = zeros(n,1);
   Z = C;
   Prob.LS.C = sparse(length(d),n);
   % Sick way to obtain the C matrix - works, but slow 
   for i=1:n
       y(i) = 1;
       Prob.LS.C(:,i) = feval(JacobMult, Z, y, -1, Prob.varargin{:});
       y(i) = 0;
   end
else
   n   = size(C,2);
   Prob.LS.C   = C;
end

m   = size(A,1);
meq = size(Aeq,1);

if length(b)~=m
   fprintf('A is %d %d. Right hand side b %d\n',size(A),length(b));
   error('A and b have not compatible dimensions');
end
if length(beq)~=meq
   fprintf('Aeq is %d %d. Right hand side beq %d\n',size(Aeq),length(beq));
   error('Aeq and beq have not compatible dimensions');
end
if m > 0 & size(A,2)~=n
   fprintf('The inequality matrix A has wrong 2nd dimension\n');
   fprintf('A is %d %d. Number of variables %d\n',size(A),n);
   return
end
if meq > 0 & size(Aeq,2)~=n
   fprintf('The equality matrix Aeq has wrong 2nd dimension\n');
   fprintf('Aeq is %d %d. Number of variables %d\n',size(Aeq),n);
   return
end

solvType=checkType('lls');

Prob.A    = [Aeq;A];
Prob.b_L  = [beq(:);-Inf*ones(m,1)];
Prob.b_U  = [beq(:);b(:)];
Prob.LS.y = d(:);
Prob.x_L  = x_L(:);
Prob.x_U  = x_U(:);
Prob.x_0  = x_0(:);
Prob.N    = n;
Prob.mLin = size(Prob.A,1);
Prob.mNonLin = 0;
Prob = tomFiles(Prob, 'ls_f', 'ls_g', 'lls_H', [], [], [], 'lls_r', 'lls_J');
Prob=ProbCheck(Prob,'lsqlin',solvType);

if isempty(Prob.SolverQP)
   Prob.SolverQP = GetSolver('lls',n > 500,1);
end

Solver=Prob.SolverQP;

% Must set SolverQP empty, if Solver is general (i.e. clsSolve), and
% might call a QP solver

Prob.SolverQP=[];
% Set default options

% TolX
Prob.optParam.eps_x=1E-6;

% TolFun
Prob.optParam.eps_f=1E-6;

% Diagnostics - off
Diagnostic=0;

% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

% LargeScale on
Prob.LargeScale=1;

% MaxIter 
Prob.optParam.MaxIter=max(400,100*n);

% PrecondBandWidth 0    TypicalX ones(numberOfVariables,1) 
% MaxPCGIter            numberOfVariables        TolPCG 0.1

if ~isempty(options)
   Prob.optParam.MaxIter = DefPar(options,'MaxIter',Prob.optParam.MaxIter);

   Prob.optParam.eps_f = DefPar(options,'TolFun',Prob.optParam.eps_f);

   if isfield(options,'LargeScale')
      if strcmpi('on',options.LargeScale)
         Prob.LargeScale=1;
      elseif strcmpi('off',options.LargeScale)
         Prob.LargeScale=0;
      end
   end
   if isfield(options,'Diagnostics')
      if strcmpi('on',options.Diagnostics)
         Diagnostic=1;
      elseif strcmpi('off',options.Diagnostics)
         Diagnostic=0;
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
   if isfield(options,'LineSearchType')
      if strcmpi('cubicpoly',options.LineSearchType)
         Prob.LineParam.LineAlg=1;
      elseif strcmpi('quadcubic',options.LineSearchType)
         Prob.LineParam.LineAlg=0;
      end
   end
end

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
           sum(~isinf(x_L)));
   fprintf(' Number of upper bound constraints:              %d\n',...
           sum(~isinf(x_U)));

   mEQ=sum(~isinf(Prob.b_U) & Prob.b_L==Prob.b_U);
   mU =sum(~isinf(Prob.b_U))-mEQ;
   mL =sum(~isinf(Prob.b_L))-mEQ;

   fprintf(' Number of linear equality constraints:          %d\n',mEQ);
   if mL > 0
      fprintf(' Number of linear lower bound constraints:       %d\n',mL);
   end
   fprintf(' Number of linear upper bound constraints:       %d\n',mU);
   fprintf('\n');
   fprintf('Algorithm selected is %s Alg # %d\n',Solver,Prob.Solver.Alg);
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf(' End diagnostic information\n');
   fprintf('\n');
end

if strcmpi('clsSolve',Solver)
   Prob.Solver.Alg=4; % Use Gauss-Newton in clsSolve
   Result = clsSolve(Prob);
else
   Result = tomRun(Solver,Prob);
end

x   = Result.x_k;
v_k = Result.v_k;
f_k = 2*Result.f_k;
r_k = Result.r_k;

Output.algorithm   = [Solver ': ' Result.SolverAlgorithm];
Output.iterations  = Result.Iter;

if nargout > 5
   Lambda.lower = zeros(n,1);
   Lambda.upper = zeros(n,1);
   ix=find(abs(x-x_L) < 1E-12);
   Lambda.lower(ix) = v_k(ix);
   ix=find(abs(x-x_U) < 1E-12);
   Lambda.upper(ix) = v_k(ix);

   Lambda.eqlin   = v_k(n+1:n+meq);
   Lambda.ineqlin = v_k(n+meq+1:n+meq+m);
   [v_k, Zv, P] = LagMult(Result.Prob,Result);

   ix = find(P(1:n) > 0);
   Lambda.upper      = zeros(n,1);
   Lambda.upper(ix)  = v_k(ix);

   ix = find(P(1:n) == -1);
   Lambda.lower      = zeros(n,1);
   Lambda.lower(ix)  = v_k(ix);
end

ExitFlag= Result.ExitFlag;
Inform  = Result.Inform;

% Convert ExitFlag to OPTIM TB
if ExitFlag==0
   ExitFlag=1;
   if PriLev > 0
      fprintf('lsqlin (%s): ',Solver);
      fprintf('Optimization Terminated successfully\n');
   end
else
   switch Inform
      case 101 
       ExitFlag=0;
       if PriLev > 0
          fprintf('lsqlin (%s): Too many iterations\n',Solver);
       end
     otherwise
       ExitFlag=-Inform;
       if PriLev > 0
          fprintf('lsqlin (%s): Possible Error.',Solver);
          ExitText = Result.ExitText;
          if ~isempty(ExitText)
             fprintf(' Solver information for Inform = %d:\n',Inform);
             for i = 1:size(ExitText,1)
                 fprintf('%s\n',ExitText(i,:));
             end
          else
             fprintf(' See solver help: Inform = %d\n',Inform);
          end
          fprintf('\n');
       end
   end
end

if nargout > 4
   Output.cgiterations= [];
   Output.firstorderopt= [];
end

% MODIFICATION LOG
%
% 011204 hkh Major revision. Add Prob as input. Selectable solver.
% 011204 hkh Test is nonempty before setting value from options into Prob
% 020128 hkh Major revision for Tomlab v4.0. Handle option JacobMult
% 040203 hkh Print ExitText if possible error (and PriLev > 0)
% 040728 med Pragma removed
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added
% 080607 hkh Use tomRun instead of tomSolve for general LLS call
% 110722 hkh Revised comments

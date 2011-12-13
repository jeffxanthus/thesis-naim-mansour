% TOMLAB LSSOL QP/LP/(LS) Matlab Solver
%
% function Result = lssolTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj.
% PriLevOpt Print Level.
% WarmStart If true, use warm start, otherwise cold start. Use WarmDefSOL
%           to set the proper parameters.
%
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution (and slacks not used) from previous run (Warm start).
% iState    Working set (if Warm start) (nb = n+nclin+ncnln) x 1 (DENSE).
%           If length(iState) < nb, setting iState(1:nb)=0;
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
% cLamda    Lagrangian multipliers from previous run (Warm start).
% SpecsFile Name of user defined SPECS file, read BEFORE optPar() is used.
% PrintFile Name of SOL Print file. Amount and type of printing determined
%           by SPECS parameters or optPar parameters.
% SummFile  Name of SOL Summary File.
% optPar    Elements > -999 take effect.
%
% -----------------------------------------------
% How optPar is used for setting SPECS parameters:
% -----------------------------------------------
%
% optPar  Structure with optimization parameters.
%
% LSSOL keywords in optPar(#):
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        10                {0,1,5,10,20,30}
%
% --- Convergence Tolerances
% 10. OPTIMALITY TOLERANCE          >0       1.1E-8            sqrt(eps)
% 11. FEASIBILITY TOLERANCE         >0       1.1E-8            sqrt(eps)
%
% --- Other Tolerances
% 21. CRASH TOLERANCE               >0       0.01      <1
% 27. RANK TOLERANCE                >0       1.1E-14           100*eps
% 30. ITERATION LIMIT               >0       max(50,5(n+m))
% 36. FEASIBILITY PHASE ITERATIONS  >0       max(50,5(n+m))
% 45. INFINITE STEP SIZE            >0       1E20
% 46. INFINITE BOUND SIZE           >0       1E20
%
% -----------------------------------------------------------------------
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from qpopt.m (similar to TOMLAB).
% Inform   LSSOL information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (LSSOL).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% xs       Solution and slack variables.
% cLamda   Lagrangian multipliers.
% iState   State for variables and constraints in iState.
%
% -----------------------------------------------------------------------
%
% For a problem description, see llsAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 2, 2000.    Last modified July 22, 2011.

function Result = lssolTL(Prob)

%#function ls_f ls_g ls_H lls_r lls_J lls_H

if nargin < 1, error('lssolTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 5; % Linear Least squares (and QP) solver

Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='LSSOL';
Result.SolverAlgorithm='LSSOL 1.05-4 QP/LP/LS code';

PriLev=Prob.PriLevOpt;

% Define lower and upper bound arrays for LSSOL
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints

BIG=1E20;

[bl, bu, n, m] = defblbu(Prob, BIG, 1);

% Initial checks on the inputs
nb = n+m;

% Check on Warm start

if Prob.WarmStart
   % Warm start for SOL solver
   Warm   = 1;
   x_0    = Prob.SOL.xs(1:n);
   if isempty(Prob.SOL.iState)
      % Use hs field
      iState = Prob.SOL.hs(1:nb);
   else
      iState = Prob.SOL.iState;
   end
   % NOT USED cLamda=Prob.SOL.cLamda;
   kx=Prob.SOL.cLamda;
else
   % Initial values
   x_0     = Prob.x_0(:);
   Warm    = 0;
   kx      = [];
   iState  = [];
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
end

probType = Prob.probType;
if checkType('lls',probType)
   % Linear least squares
   H = full(Prob.LS.C);
   y = full(Prob.LS.y);
   c = [];
elseif checkType('lp',probType) | checkType('qp',probType)
   % Check if any linear part
   c = full(Prob.QP.c(:));
   H = full(Prob.QP.F);
   y = [];
else
   fprintf('\n\nCalling LSSOL with ');
   fprintf('problem type %d: %s\n\n',probType,checkType([],probType));
   error('LSSOL cannot handle this problem type, only handles lp,qp,lls')
end

nnObj = size(H,1); % number of nonlinear variables

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   if mA~=m, error('Linear constraints A MUST have m rows!'); end
end

% Determine type of problem
if ~isempty(y)
   % Linear least squares
   r = H*x_0-y;
   Result.f_0=0.5*(r'*r);
elseif isempty(c) | all(c==0)
   if isempty(H) | nnz(H) == 0
      Result.f_0=0;
   else
      Result.f_0=0.5*(x_0(1:nnObj)'*H*x_0(1:nnObj));
   end
else
   if isempty(H) | nnz(H) == 0
      Result.f_0=c(1:n)'*x_0(1:n);
   else
      Result.f_0=0.5*(x_0(1:nnObj)'*H*x_0(1:nnObj)) + c(1:n)'*x_0(1:n);
   end
end

% Set up the constraint matrix A
if m==0     % isempty(Prob.A)
   Prob.A = zeros(0,n);
end

if PriLev > 2
   fprintf('# of non zero entries in linear constraints %d. ',nnz(Prob.A))
   if PriLev > 5
      fprintf('Number of variables          %d.\n',n)
      fprintf('Number of linear constraints %d.\n',m)
   end
   if PriLev >= 1000
      xprinte(bl(1:n),'x_L: ');
      xprinte(bu(1:n),'x_U: ');
      xprinte(bl(n+1:nb),'b_L: ');
      xprinte(bu(n+1:nb),'b_U: ');
      xprinte(c,'c: ');
      if PriLev > 0, pause; end
   end
   if PriLev >= 7
      if PriLev >= 10, pause; end
      PrintMatrix(full(Prob.A),'Constraint matrix A:')
      if PriLev >= 10, pause; end

      disp('x_L x_0 x_U');
      mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 10, pause; end
   end
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

Upper    = 0;
optPar   = Prob.SOL.optPar(:)';

[x_k, Inform, iState, cLamda, Iter, Obj, rProj, kx] = ...
   lssol(full(Prob.A), bl, bu, c, x_0, optPar, H, y, Warm, ...
         iState, Upper, kx, ...
         Prob.SOL.SpecsFile, Prob.SOL.PrintFile, Prob.SOL.SummFile, ...
         PriLev, Prob.Name);

switch Inform
   case 4
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 3
     ExitFlag=4;  % Infeasible
   case 1
     ExitFlag=3;  % Rank problem
   case 6
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k = Obj;
% Determine type of problem
if ~isempty(y)
   % Linear least squares
   Result.r_k = H*x_k-y;
end
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;
% Change sign when iState == 2
ix             = find(iState==2);
Result.v_k(ix) = -cLamda(ix);

% Warm start
% The dummy objective row is last in xs, value -Obj

% Could use Result.x_k instead
Result.SOL.xs=[x_k;zeros(m,1);-Obj];

Result.SOL.iState=iState;
% Put iState also in hs field
Result.SOL.hs=[iState;0];

% Could use Result.v_k instead
%Result.SOL.cLamda=cLamda;
% PUT kx here for warm start
Result.SOL.cLamda=kx;
%format long
%[cLamda iState]

Result.SOL.optPar=optPar;

if ~isempty(y)
   Result.g_k=H'*(H*x_k-y);
elseif ~isempty(c)
   if ~isempty(H)
      Result.g_k=H*x_k+c;
   else
      Result.g_k=c;
   end
elseif isempty(H)
   Result.g_k=[];
else
   Result.g_k=H*x_k;
end

Result.FuncEv    = Iter;
Result.GradEv    = Iter;
Result.ConstrEv  = Iter;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

%Compute Result.xState and Result.QP.B only
Result = StateDef(Result, x_k, [], [], Prob.optParam.xTol, [], [], bl, bu, 1);

% Use already computed iState
Result.bState                          = iState(n+1:n+m);

switch Inform
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'Weak local solution (nonunique) was reached';
   case 2
      Text = 'The solution appears to be unbounded';
   case 3
      Text = str2mat('The constraints could not be satisfied.' ...
                    ,'The problem has no feasible solution');
   case 4
      Text = 'Too many iterations, in either phase';
   case 5
      Text = str2mat('50 iterations of working set without change in x' ...
             ,'Cycling?');
   case 6
      Text = 'An input parameter was invalid';
   otherwise
      Text = 'NOTE: UNKNOWN LSSOL Inform value.';
end

Result.ExitText = Text;


if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLSSOL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LSSOL: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');


   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector iState for x and constraints = \n');
      xprinti(iState(1:min(length(iState),MAX_x)),'iState: ');
   end

   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (cLamda) = \n');
      xprinte(cLamda(1:min(length(cLamda),MAX_c)),'cLamda:');
   end
end

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 001102 hkh  New lssolMex interface.
% 010710 hkh  Error in comments
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 020417 hkh  Send optPar back in Result.SOL
% 020821 hkh  Remove code appearing twice
% 030508 hkh  Weight problem according to Tomlab standard using WeightType
% 030521 hkh  For LLS, compute r_k = H*x_k-y, LSSOL returns projected r_k
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040105 hkh  No longer necessary to define 0-matrix if LP problem
% 040126 hkh  Adding check if problemType is lls, lp or qp
% 040602 med  Help fix PriLev to PriLevOpt
% 041201 hkh  weightType=3, dynamic weights, not possible, stop with error
% 041202 hkh  Revise call to defblbu and StateDef
% 041222 med  Safeguard added for x_0
% 050614 med  Updated help for optPar vector
% 050614 med  Updated iState, removed optPar 3 and 4
% 050617 hkh  Removed SOLSet call
% 050617 hkh  Use iState to compute StateDef fields, B and xState
% 050725 med  Modified Prob.A no longer passed to Result
% 050829 med  Updated StateDef call
% 060818 hkh  Skip using dummy if no linear constraint, set A=zeros(0,n)
% 060818 med  isnan checks removed for x_0
% 070621 frhe all(H == 0) changed to nnz(H) == 0 to avoid vector output
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
% 090717 med  Fixed f calculations
% 090810 hkh  Avoid a copy A of Prob.A
% 110722 hkh  Removed weightType and weightY code, weighting done in llsAssign

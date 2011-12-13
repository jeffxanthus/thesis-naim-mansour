% TOMLAB LPOPT LP Solver
%
% function Result = lpoptTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format.
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
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
% LPOPT keywords in optPar(#):
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        10                {0,1,5,10,20,30}
%
% --- Convergence Tolerances
%     In LPOPT/QPOPT: macheps = 2^(-53);  eps in Matlab is = 2^(-52); 
% 10. OPTIMALITY TOLERANCE          >0       1.0537E-8         sqrt(macheps)
% 11. FEASIBILITY TOLERANCE         >0       1.0537E-8         sqrt(macheps)
%
% --- Other Tolerances
% 21. CRASH TOLERANCE               >0       0.01      <1
% 27. RANK TOLERANCE                >0       1.1102E-14        100*macheps
% 30. ITERATION LIMIT               >0       max(2000,5(n+m))
% 33. MIN SUM YES (or NO)           0        0         1       1=min infeas
%     If 1 (MIN SUM YES), minimize the infeasibilities before return
% 36. FEASIBILITY PHASE ITERATIONS  >0       max(2000,5(n+m))
% 45. INFINITE STEP SIZE            >0       1E20
%
% --- Frequencies
% 51. CHECK FREQUENCY               >0       50
% 52. EXPAND FREQUENCY              >0       5
%
% -----------------------------------------------------------------------
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient (trivially c).
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from lpopt.m (similar to TOMLAB).
% Inform   LPOPT information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (LPOPT).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% xs       Solution and slack variables.
% iState   State for variables and constraints in iState.
%
% -----------------------------------------------------------------------
%
% For a problem description, see lpAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 25, 2000.  Last modified Jun 6, 2008.

function Result = lpoptTL(Prob)

if nargin < 1, error('lpoptTL needs the Prob structure as input'); end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 1; % LP solver

Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='LPOPT';
Result.SolverAlgorithm='LPOPT 1.0-10 LP code';

BIG=1E20;

[bl, bu, n, m] = defblbu(Prob, BIG, 1);

PriLev=Prob.PriLevOpt;

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
      iState=Prob.SOL.iState;
   end
   % NOT USED cLamda=Prob.SOL.cLamda;
else
   Warm   = 0;
   x_0    = Prob.x_0(:);
   iState = [];
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
end

% Check if any linear part
c = Prob.QP.c(:);

% Determine type of problem
if isempty(c) | all(c==0)
   Result.f_0=0;
else
   Result.f_0=c(1:n)'*x_0(1:n);
end

% Set up the constraint matrix A 
if m==0     % isempty(Prob.A)
   A=zeros(0,n);
else
   A = Prob.A;
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

%[optPar,SpecsFile,PrintFile,SummFile]=SOLSet('lpopt',8,0,0,size(Prob.A,1),Prob);

% Defaults in the LPOPT and QPOPT code
% Currently no difference between LP and QP problems
% macheps in QPOPT is different from SNOPT, MINOS and Matlab
% macheps = 2^(-53);  % eps in Matlab is = 2^(-52); 
% tolfea = sqrt(macheps);  % 1.05E-8, epspt5
% tolOpt = macheps^(-0.8); % 1.72E-13, epspt8
% tolrnk = 100*macheps;    % 1.11E-14

% 1.  PRINT LEVEL  optPar(1)  = 10;
% 10. OPTIMALITY TOLERANCE   optPar(10) = 1.72E-13; (macheps^(-0.8)), epspt8
% 11. FEASIBILITY TOLERANCE  optPar(11) = 1.05E-8; (sqrt(macheps))
% 21. CRASH TOLERANCE        optPar(21) = 0.01; 
% 27. RANK TOLERANCE         optPar(27) = 1.11E-14; (100*macheps)
% 30. ITERATION LIMIT        optPar(30) = max(50,5*(n+m));
% 33. MIN SUM. IF SET TO 1, minimize infeasibilities before return (def 0)
% 36. FEAS PHASE ITERATIONS LIMIT  optPar(36) = max(50,5*(n+m));
% 45. INFINITE STEP SIZE optPar(45) = 1E20;
% 48. MAX DEGREES OF FREEDOM, ONLY USED IF HESSIAN ROWS == N (Def = n)
% 51. CHECK FREQUENCY  optPar(51) = 50;
% 52. EXPAND FREQUENCY optPar(52) = 5;
% 3,4,47 set by Matlab interface
% 3.  PRINT FILE           0        0                 Fortran Unit #
%           SET BY INTERFACE IF PrintFile is given
% 4.  SUMMARY FILE         0        0                 Fortran Unit #
%           SET BY INTERFACE IF SummFile  is given
% 47. HESSIAN ROWS         0        n         n       0 if FP or LP
%     IMPLICITLY GIVEN BY THE DIMENSIONS OF H IN THE CALL FROM MATLAB

optPar                      = Prob.SOL.optPar(:)';
% Length 52 used in MEX, but defined optPar length with 62 elements in MEX
optPar(length(optPar)+1:52) = -900;

if optPar(10) <= 0
   optPar(10) = sqrt(2^(-53));
end
if optPar(30) <= 0
   optPar(30) = max(2000, 5*nb); 
end
if optPar(36) <= 0
   optPar(36) = max(2000, 5*nb); 
end

[Inform, Iter, iState, Ax, cLamda, Obj, x_k] = lpopt( ...
         full(A), bl, bu, c, Warm, x_0, iState, ...
         Prob.SOL.SpecsFile, Prob.SOL.PrintFile, Prob.SOL.SummFile, ...
         PriLev, optPar );

switch Inform
   case 4
     ExitFlag=1;  % Too many iterations
   case 2
     ExitFlag=2;  % Unbounded
   case 3
     ExitFlag=4;  % Infeasible
   case {1}
     ExitFlag=3;  % ??? May be optimal
   case {5,6,7}
     ExitFlag=10; % Input errors
   otherwise
     ExitFlag=0;
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;
Result.g_k = c;

% Warm start
% The dummy objective row is last in xs, value -Obj

% Could use Result.x_k instead
Result.SOL.xs=[x_k;zeros(m,1);-Obj];

Result.SOL.iState=iState(1:nb);
% Put iState also in hs field
Result.SOL.hs=[iState;0];

% Could use Result.v_k instead
%Result.SOL.cLamda=cLamda;

Result.SOL.optPar=optPar;

Result.FuncEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState and Result.QP.B only
Result = StateDef(Result, x_k, [], [], Prob.optParam.xTol, [], [], bl, bu, 1); 

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%Result.xState=iState(1:n);
% Use already computed iState
Result.bState=iState(n+1:n+m);

switch Inform
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'A dead point was reached';
   case 2
      Text = 'The solution appears to be unbounded (or badly scaled)';
   case 3
      Text = str2mat('The constraints could not be satisfied.' ...
                    ,'The problem has no feasible solution');
   case 4
      Text = 'Too many iterations, in either phase';
   case 5
      Text = str2mat('The Maximum degrees of freedom is too small' ...
             ,'The reduced Hessian must expand if further progress is' ...
             ,'too be made');
   case 6
      Text = 'An input parameter was invalid';
   case 7
      Text = 'The problem type was not recognized';
   otherwise
      Text = 'NOTE: UNKNOWN LPOPT Inform value.';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLPOPT solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LPOPT: Inform = %2d, ',Inform)
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
% 000716 hkh New lpoptMex interface.
% 000916 hkh Return string ExitText with interpretation of Inform flag 
% 001004 hkh Revision for direct call to the mex file
% 010710 hkh Error in comments
% 020417 hkh Send optPar back in Result.SOL
% 020821 hkh Remove code appearing twice
% 040102 hkh Revision for v4.2, call iniSolve and endSolve
% 040102 hkh Return only Iter ~= 0, FuncEv=ConstrEv=0
% 040602 med Help fix PriLev to PriLevOpt
% 041202 hkh Revise calls to defblbu and StateDef
% 041222 med Safeguard added for x_0
% 050614 med Updated help for optPar vector
% 050614 med Updated iState, removed optPar 3 and 4
% 050725 med Modified Prob.A no longer passed to Result
% 060818 hkh Skip using dummy if no linear constraint, set A=zeros(0,n)
% 060818 med isnan checks removed for x_0
% 070524 med Removed print out, return statement and size(A)
% 080213 hkh Removed SOLSet, set only optPar(10,30,36)
% 080606 med Switched to iniSolveMini
% 080607 hkh Switched to endSolveMini

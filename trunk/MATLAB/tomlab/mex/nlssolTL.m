% TOMLAB NLSSOL Constrained nonlinear least squares solver
%
% function Result = nlssolTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format.
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% c_L, c_U  Bounds on nonlinear constraints.
% A         Linear constraint matrix.
% PriLevOpt Print level.
% WarmStart If true, use warm start, otherwise cold start. Use WarmDefSOL
%           to set the proper parameters.
%
% -----------------------------------------------
% Fields used in Prob.SOL:
% -----------------------------------------------
% xs        Solution from previous run, elements xs[1:n].
% iState    Working set (if Warm start) (nb = n+nclin+ncnln) x 1 (DENSE).
%           If length(iState) < nb, setting iState(1:nb)=0;
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
% cLamda    Lagrangian multipliers (dual solution vector) (nb x 1 vector).
% H         Cholesky factor of Hessian approximation.
%           Hessian no  - reordered variables.
%           Hessian yes - natural order of variables, used for Warm start.
%
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
% NLSSOL keywords in optPar(#):
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by NLSSOL,
% if not the value is altered in the SPECS file (input SpecsFile)
%
% # is the index in the optPar vector.
% See NPSOL User's Guide for the input (SPECS) keywords and description
% Additional keyword for NLSSOL are
% JTJ Initial Hessian (default) or Unit Initial Hessian
% Reset frequency (default 2)
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        10                {0,1,5,10,20,30}
%     or MAJOR PRINT LEVEL
% 2.  MINOR PRINT LEVEL             0        0                 {0,1,5,10,20,30}
%
% --- Convergence Tolerances
% 9.  NONLINEAR FEASIBILITY TOLERANCE >0     1.1E-8            sqrt(eps)
% 10. OPTIMALITY TOLERANCE          >0       3.0E-13           eps^0.8
% 11. LINEAR FEASIBILITY TOLERANCE  >0       1.1E-8            sqrt(eps)
%
% --- Derivative checking
% 13. VERIFY LEVEL                  -1       -1        3       {-1,0,1,2,3}
% 14. START OBJECTIVE CHECK AT COL  0        1         n
% 15. STOP OBJECTIVE CHECK AT COL   0        n         n
% 16. START CONSTRAINT CHECK AT COL 0        1         n
% 17. STOP CONSTRAINT CHECK AT COL  0        n         n
%
% --- Other Tolerances
% 21. CRASH TOLERANCE               >0       0.01      <1
% 22. LINESEARCH TOLERANCE          >0       0.9       <1
%     Line search accuracy
%
% 30. ITERATION LIMIT               >0       max(50,3(n+m_L)+10*m_N)
% 36. MINOR ITERATIONS LIMIT        >0       max(50,3(n+m_L+m_N))
% 37. STEP LIMIT                    >0       2
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     Is set by nlssolTL dependent on Prob.ConsDiff, Prob.NumDiff
% 41. FUNCTION PRECISION            >0       3.0E-13           eps^0.8=eps_R
%     Relative tolerance on function
% 42. DIFFERENCE INTERVAL           >0       5.48E-8           eps^0.4
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.70E-5           eps^{0.8/3}
% 45. INFINITE STEP SIZE            >0       max(BIGBND,1E10)
% 46. INFINTE BOUND SIZE            >0       1E10              = BIGBND
% 47. JTJ INITIAL HESSIAN            0       1         1       0 = UNIT
% 48. RESET FREQUENCY                0       2                 Reset of (47)
% 50. HESSIAN YES or NO              0       0                 1 = YES
%
% -----------------------------------------------------------------------
%
% OUTPUT:
%
% Result     Structure with results (see ResultDef.m):
%
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Gradient of the function.
% r_k      Residual vector.
% J_k      Jacobian matrix.
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% rc       Reduced costs. If ninf=0, last m == -v_k.
% ExitFlag Exit status from nlssol.m (similar to TOMLAB).
% Inform   NLSSOL information parameter.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (nlssol).
% SolverAlgorithm  Description of the solver.
%
% The following outputs are set in the Result.SOL sub field:
%
% x        Solution vector (n by 1) with n decision variable values
% iState   Basis status of constraints + variables, (m + n x 1 vector)
%          State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%          2=superbasic (between bounds), 3=basic (between bounds).
% cLamda   Lagrangian multipliers (dual solution vector) (m x 1 vector).
% H        Cholesky factor of Hessian approximation.
%          Hessian no  - reordered variables.
%          Hessian yes - natural order of variables, used for Warm start.
%
% -----------------------------------------------------------------------
%
% For a problem description, see clsAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Sept 17, 2000.  Last modified Jul 17, 2009.

function Result = nlssolTL(Prob)

%#function ls_f ls_g ls_H lls_r lls_J lls_H ls_rJ

if nargin < 1, error('nlssolTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

if isempty(Prob.LS)
   fprintf('Prob.LS not defined\n');
   error('No least squares problem')
end
%if any(Prob.probType ==[2 7 8 11])
%   error('nlssol cannot solve LP, MILP, QP and MIQP problems');
%end

Prob.solvType = 6; % Constrained NLLS solver

if isempty(Prob.NumDiff)
   ObjDers = 1;
elseif Prob.NumDiff==6
   ObjDers = 0;
else
   ObjDers = 1;
end
if Prob.mNonLin > 0
   if isempty(Prob.ConsDiff)
      ConsDers = 1;
   elseif Prob.ConsDiff==6
      ConsDers = 0;
   else
      ConsDers = 1;
   end
else
   ConsDers = 0;
end

Prob = iniSolve(Prob,6,ObjDers,ConsDers);

Result=ResultDef(Prob);
Result.Solver='NLSSOL';
Result.SolverAlgorithm='NLSSOL 5.02 constrained NLLS code';

PriLev=Prob.PriLevOpt;

% Define lower and upper bound arrays for NLSSOL
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints

BIG=1E20;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

m  = m1+m2;
nb = n+m1+m2;

% Initial checks on the inputs

% Check on Warm start, then set iState,x,H and cLamda

if Prob.WarmStart
   % Warm start for dense SOL solver
   Warm = 1;
   x_0    = Prob.SOL.xs(1:n);
   iState = Prob.SOL.iState;
   if length(iState) < nb
      % Use hs field
      iState = Prob.SOL.hs(1:nb);
   end
   cLamda = Prob.SOL.cLamda;
   H      = Prob.SOL.H;
else
   % Initial values
   Warm = 0;
   x_0    = Prob.x_0(:);
   if length(x_0) < n, x_0=zeros(n,1); end
   x_0    = max(bl(1:n),min(x_0,bu(1:n)));
   iState = []; H = []; cLamda = [];
end

nEqual = 0;

if m2 > 0
   % How many nonlinear equality constraints
   nEqual=sum(bl(n+m1+1:n+m)==bu(n+m1+1:n+m));

   % Must call function, in case global values are needed for constraints
   % Because nlssol calls constraints first, then functions
   r_0        = nlp_r(x_0, Prob);
   Result.r_0 = r_0;
   Result.f_0 = 0.5*(r_0'*r_0);
   % Sometimes maybe also a call to the gradient might be needed.
   % Do not make such a call now, it is highly unlikely.
else
   r_0        = nlp_r(x_0, Prob);
   Result.r_0 = r_0;
   Result.f_0 = 0.5*(r_0'*r_0);
   %Result.f_0 = c'*x_0;
end

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end
end

if PriLev > 2
   if PriLev >= 3
      PrintMatrix(full(Prob.A),'Constraint matrix A:')
      if PriLev >= 10, pause; end
   end

   if PriLev >= 1000
      xprinte(bl(1:n),'x_L: ');
      xprinte(bu(1:n),'x_U: ');
      xprinte(bl(n+1:nb),'b_L: ');
      xprinte(bu(n+1:nb),'b_U: ');
      if PriLev > 0, pause; end
   end
   if PriLev > 5
      fprintf('Non linear variables %d. ',n)
      if m2 > 0
         fprintf('Non linear constraints %d. ',m2)
      end
      fprintf('Equalities %d.\n',nEqual)
      fprintf('Total number of constraints %d.\n',m1+m2)
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

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('nlssol',6, n,m2,m1,Prob);

%     % NLSSOL. Difference NPSOL is #47 and #48
%     % nnObj == n, nnJac == nonlin cons (ncnln), m = linear cons (nclin)
%
%     % Currently no difference between LLS, LP and QP problems
%
%     seps=sqrt(eps);                    % 1.1E-8
%     e08 = 3.0E-13;                     % eps^(0.8)
%     m30=max(50,3*(nObj+m) + 10*nJac);  % 3*(n+nclin) + 10*ncnln
%     m36=max(50,3*(nObj + m + nJac));   % 3*(n+nclin+ncnln)
%     e09 = eps^(0.9);
%     %      [   1     2    3    4    5    6    7    8     9    10
%     optPar=[   0     0   -1   -1 -999 -999 -999 -999  seps   e08 ...
%             seps  -999    0    1 nObj    1 nObj -999  -999  -999 ...
%             0.01   0.9 -999 -999 -999 -999 -999 -999  -999   m30 ...
%             -999  -999 -999 -999 -999  m36 -999 -999     3  -999 ...
%              e08  -900 -900 -999 1E10 1E10    1    2  -999  -999 ...
%             -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
%             -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
%             -999  ...
%            ];

optPar    = Prob.SOL.optPar(:)';
optPar(length(optPar)+1:50)=-999;

if Warm > 0
   % Warm start for dense SOL solver, must set Hessian Yes !!!
   optPar(50)=1;
else
   %optPar(50)=0;
end

if Prob.NumDiff == 6
   if Prob.ConsDiff == 6
      optPar(39)=0;
   else
      optPar(39)=2;
   end
elseif Prob.ConsDiff == 6
   optPar(39)=1;
else
   optPar(39)=3;
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

% y = Prob.LS.y;
% Use clean formulation, without y
m = size(Prob.LS.y,1);

% Safe guard
if m == 1 & size(Prob.LS.y,2) > 1
   m = size(Prob.LS.y,2);
end
y = zeros(m,1);
% Avoid derivative check
if optPar(13) < 0
   optPar(13)=-1;
end

[x, Inform, iState, cLamda, iwCount, fObj, gObj, r, J, fCon, gCon, H] = ...
    nlssol( full(Prob.A), bl, bu, x_0, Prob, y, optPar, ...
           Warm, H, iState, cLamda, ...
           Prob.SOL.SpecsFile,Prob.SOL.PrintFile,Prob.SOL.SummFile,...
           PriLev, Prob.Name );

switch Inform
   case 4
     ExitFlag=1;  % Too many iterations
   %case ?
   %  ExitFlag=2;  % Unbounded
   case {2,3}
     ExitFlag=4;  % Infeasible
   case {1,6}
     ExitFlag=3;  % Rank problem
   case {7,9}
     ExitFlag=10; % Input errors
   case 0
     ExitFlag=0;  % Optimal solution found
   otherwise
     ExitFlag=10;
end

% Result.r_k = y-r;
Result.r_k = nlp_r(x,Prob);
Result.J_k = nlp_J(x,Prob);
Result.f_k = fObj;
Result.g_k = gObj;
Result.c_k = fCon;
Result.x_k = x;
Result.x_0 = x_0;
Result.v_k = cLamda;

if ~isempty(gCon)
   % TOMLAB now has the format one constraint / row, same as SOL solvers
   Result.cJac=gCon;
end

% Saved for a warm start
% The dummy objective row is last in xs, value -fObj

% Could use Result.x_k instead
Result.SOL.xs=[x;zeros(m,1);-fObj];

Result.SOL.iState=iState;
% Put iState also in hs field
Result.SOL.hs=[iState;0];

Result.SOL.H=H;

% Could use Result.v_k instead
Result.SOL.cLamda=cLamda;

Result.SOL.optPar=optPar;

Result.Iter      = iwCount(1);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

optParam = Prob.optParam;
if m1 > 0
   Result = StateDef(Result, x(1:n), Prob.A*x, fCon, ...
                  optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);
else
   Result = StateDef(Result, x(1:n), [], fCon, ...
                  optParam.xTol, optParam.bTol, optParam.cTol, bl, bu, 1);
end

switch Inform
   case 0
      Text = 'Optimal solution found';
   case 1
      Text = str2mat('Optimal solution found' ...
                    ,'but not to requested accuracy');
   case 2
      Text = 'No feasible point for the linear constraints';
   case 3
      Text = 'No feasible point for the nonlinear constraints';
   case 4
      Text = 'Too many major iterations';
   case 5
      Text = 'Problem is unbounded, or badly scaled';
   case 6
      Text = 'The current point cannot be improved on';
   case 7
      Text = 'Large errors found in the derivatives';
   case 9
      Text = str2mat('Errors found in the input parameters ' ...
             ,'Problem abandoned');
   otherwise
      Text = 'User requested termination';
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nNLSSOL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('NLSSOL: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',fObj);
   fprintf('Iterations%7d. ',iwCount(1));
   fprintf('\n');

   fprintf('fObj and gObj evaluations%7d %d\n',iwCount(2:3));

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=length(x);
      end
      fprintf('Optimal x = \n');
      xprinte(x(1:min(n,MAX_x)),'x:  ');
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

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 000917 hkh  New nlssolMex interface.
% 011212 hkh  Use Prob.Name instead of creating new variable ProbName
% 020304 hkh  Set optPar(39) DERLVL dependent on Prob.ConsDiff, Prob.NumDiff
% 020417 hkh  Send optPar back in Result.SOL
% 030902 ango Updated comments
% 031126 hkh  Call nlp_r and compute objective instead of calling nlp_f
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040524 hkh  Use nlp_J to get Result.J_k, nlssol sometimes returns non-optimal
% 040524 hkh  Use nlp_r to get Result.r_k, seems nlssol sometimes returns -r
% 040602 med  Help fix PriLev to PriLevOpt
% 041202 hkh  If Prob.LS not defined, stop execution with error message
% 041202 hkh  Revise call to defblbu and StateDef
% 041202 hkh  If internal differentiation, different call to iniSolve
% 041222 med  Safeguard added for x_0s
% 050306 hkh  Step limit added as control parameter
% 050614 med  Updated help for optPar vector
% 050614 med  Updated iState, removed optPar 3 and 4
% 051216 med  optPar(13) default to -1 (no derivative check)
% 080276 hkh  Avoid call to SOLSet
% 080604 hkh  Safeguard call to iniSolve for empty NumDiff and ConsDiff
% 090717 med  Residual calculation updated

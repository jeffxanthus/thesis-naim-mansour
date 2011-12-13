% TOMLAB NPSOL NLP Solver
%
% function Result = npsolTL(Prob)
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
% NPSOL keywords in optPar(#):
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by NPSOL,
% if not the value is altered in the SPECS file (input SpecsFile)
%
% # is the index in the optPar vector.
% See NPSOL User's Guide for the input (SPECS) keywords and description
%
% #   SPECS keyword text            Lower    Default   Upper   Comment
%
% --- Printing
% 1.  PRINT LEVEL                   0        10                {0,1,5,10,20,30}
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
%                            Default 1E-10 if any pair of finite bounds, n<1000
% 22. LINESEARCH TOLERANCE          >0       0.9       <1
% 30. ITERATION LIMIT               >0       max(50,3(n+m_L)+10*m_N)
% 36. MINOR ITERATIONS LIMIT        >0       max(50,3(n+m_L+m_N))
% 37. STEP LIMIT                    >0       2
% 39. DERIVATIVE LEVEL              0        3         3       {0,1,2,3}
%     Is set by npsolTL dependent on Prob.ConsDiff, Prob.NumDiff
% 41. FUNCTION PRECISION            >0       3.0E-13           eps^0.8=eps_R
% 42. DIFFERENCE INTERVAL           >0       5.48E-8           eps^0.4
% 43. CENTRAL DIFFERENCE INTERVAL   >0       6.70E-5           eps^{0.8/3}
% 45. INFINITE STEP SIZE            >0       max(BIGBND,1E10)
% 46. INFINITE BOUND SIZE           >0       1E10              = BIGBND
% 50. HESSIAN YES or NO              0       0                 1 = YES
%
% -----------------------------------------------------------------------
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Gradient of the function.
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from npsol.m (similar to TOMLAB).
% Inform   NPSOL information parameter.
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (npsol).
% SolverAlgorithm  Description of the solver.
%
% The following output are set in the Result.SOL sub field
% Warm start
% x         Solution vector (n by 1) with n decision variable values.
% iState    Working set (if Warm start) (nb = n+nclin+ncnln) x 1 (DENSE).
%           If length(iState) < nb, setting iState(1:nb)=0;
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
% cLamda    Lagrangian multipliers (dual solution vector) (m x 1 vector).
% H         Cholesky factor of Hessian approximation.
%           Hessian no  - reordered variables.
%           Hessian yes - natural order of variables, used for Warm start.
%
% -----------------------------------------------------------------------
%
% For a problem description, see conAssign.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sept 16, 2000.  Last modified Oct 2, 2009.

function Result = npsolTL(Prob)

%#function nlp_cdc nlp_fg

if nargin < 1, error('npsolTL needs the Prob structure as input');end

% funfdf Name of routine [f,gradf] = funfdf(x, Prob, mode, nstate)
%        funfdf=nlp_fg, included in TOMLAB.
% funcdc Name of routine [g,gJac]  = funcdc(x, Prob, mode, nstate, needg)
%        funcdc=nlp_cdc, included in TOMLAB.

global MAX_x MAX_c% Max number of variables/constraints to print

Prob.solvType = 3; % NLP (CON) solver

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

Prob = iniSolve(Prob,3,ObjDers,ConsDers);

Result=ResultDef(Prob);
Result.Solver='NPSOL';
Result.SolverAlgorithm='NPSOL 5.02 NLP code';

PriLev=Prob.PriLevOpt;

% Define lower and upper bound arrays for NPSOL
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

% Initial checks on the inputs
m       = m1+m2;

nb     = n+m1+m2;

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
   % Because npsol calls constraints first, then functions
   Result.f_0= nlp_f(x_0, Prob);
   % Sometimes maybe also a call to the gradient might be needed.
   % Do not make such a call now, it is highly unlikely.
else
   Result.f_0= nlp_f(x_0, Prob);
   %Result.f_0 = c'*x_0;
end

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end 
end 

if PriLev > 2 
   if PriLev >= 3 & ~isempty(Prob.A)
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
      if PriLev >= 1000, pause; end
      PrintMatrix(full(Prob.A),'Constraint matrix A:')
      if PriLev >= 1000, pause; end
   
      disp('x_L x_0 x_U');
      %xprint(x_0,'x_0: ',' %14.9f',5);
      %xprinte(bl,'bl: ');
      %xprinte(bu,'bu: ');
      mPrint([bl(1:n) x_0 bu(1:n)],'xl/x0/xu:')
      mPrint([bl(n+1:n+m) bu(n+1:n+m)],'bl/bu:')
      if PriLev >= 1000, pause; end
   end
end

%[optPar, SpecsFile, PrintFile, SummFile] = SOLSet('npsol',3,n,m2,m1,Prob);

%      % nnObj == n, nnJac == nonlin cons (ncnln), m = linear cons (nclin)
%
%      seps=sqrt(eps);                    % 1.1E-8
%      e08 = 3.0E-13;                     % eps^(0.8)
%      m30=max(50,3*(nObj+m) + 10*nJac);  % 3*(n+nclin) + 10*ncnln 
%      m36=max(50,3*(nObj + m + nJac));   % 3*(n+nclin+ncnln) 
%      e09 = eps^(0.9);
%      %      [   1     2    3    4    5    6    7    8     9    10
%      optPar=[   0     0   -1   -1 -999 -999 -999 -999  seps   e08 ...
%              seps  -999    0    1 nObj    1 nObj -999  -999  -999 ...
%              0.01   0.9 -999 -999 -999 -999 -999 -999  -999   m30 ...
%              -999  -999 -999 -999 -999  m36 -999 -999     3  -999 ...
%               e08  -900 -900 -999 1E10 1E10 -999 -999  -999     0 ...
%              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
%              -999  -999 -999 -999 -999 -999 -999 -999  -999  -999 ...
%              -999  ...
%             ];

optPar    = Prob.SOL.optPar(:)';
optPar(length(optPar)+1:50)=-999;
if n < 1000 & optPar(21)==-999
   if any(isfinite(bl(1:n))&isfinite(bu(1:n)))
      optPar(21)=1E-10;
   end
end

if Warm > 0 
   % Warm start for dense SOL solver, must set Hessian Yes !!!
   optPar(50)=1;
else
   %optPar(50)=0;
end

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
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
% Avoid derivative check
if optPar(13) < 0
   optPar(13)=-1;
end

[x, Inform, iState, cLamda, iwCount, fObj, gObj, fCon, gCon, H] = ...
    npsol( full(Prob.A), bl, bu, x_0, Prob, optPar, ...
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
   fprintf('\nNPSOL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('NPSOL: Inform = %2d, ',Inform)
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
if any(Inform == [0;1;6])
   Result=cErrCompute(Result);
end

% Check of correct return codes from NPSOL
function Result = cErrCompute(Result)

% 9.  NONLINEAR FEASIBILITY TOLERANCE >0     1.1E-8            sqrt(eps)
% 10. OPTIMALITY TOLERANCE          >0       3.0E-13           eps^0.8
% 11. LINEAR FEASIBILITY TOLERANCE  >0       1.1E-8            sqrt(eps)

bTol = Result.Prob.SOL.optPar(11);
cTol = Result.Prob.SOL.optPar(9);
if bTol < -998.9
   bTol = 1.1e-8;
end
if cTol < -998.9
   cTol = 1.1e-8;
end

x_k  = Result.x_k;
if isempty(x_k)
   return
end
Prob = Result.Prob;
bErr = 0;
cErr = 0;
if Prob.mLin > 0
   if isempty(Result.Ax)
      Ax = Prob.A * x_k;
      Result.Ax = Ax;
   else
      Ax = Result.Ax;
   end
   ixE  = find(Prob.b_L == Prob.b_U);
   if ~isempty(ixE)
      bErr = max(bErr, max(abs(Prob.b_U(ixE)-Ax(ixE))));
   end
   ixL  = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_L));
   if ~isempty(ixL)
      bErr = max(bErr, max(max(0,Prob.b_L(ixL)-Ax(ixL))));
   end
   ixU  = find(Prob.b_L ~= Prob.b_U & ~isinf(Prob.b_U));
   if ~isempty(ixU)
      bErr = max(bErr, max(max(0,Ax(ixU)-Prob.b_U(ixU))));
   end
end
if Prob.mNonLin > 0
   c_k  = Result.c_k;
   ixE  = find(Prob.c_L == Prob.c_U);
   if ~isempty(ixE)
      cErr = max(cErr, max(abs(Prob.c_U(ixE)-c_k(ixE))));
   end
   ixL  = find(Prob.c_L ~= Prob.c_U & ~isinf(Prob.c_L));
   if ~isempty(ixL)
      cErr = max(cErr, max(max(0,Prob.c_L(ixL)-c_k(ixL))));
   end
   ixU  = find(Prob.c_L ~= Prob.c_U & ~isinf(Prob.c_U));
   if ~isempty(ixU)
      cErr = max(cErr, max(max(0,c_k(ixU)-Prob.c_U(ixU))));
   end
end

if bErr > bTol
   Result.ExitText = 'No feasible point for the linear constraints';
   Result.Inform = 2;
   Result.ExitFlag = 4;
elseif cErr > cTol
   Result.ExitText = 'No feasible point for the nonlinear constraints';
   Result.Inform = 3;
   Result.ExitFlag = 4;
end

% MODIFICATION LOG:
%
% 000916 hkh New npsolMex interface.
% 011212 hkh Use Prob.Name instead of creating new variable ProbName
% 020304 hkh Set optPar(39) DERLVL dependent on Prob.ConsDiff, Prob.NumDiff
% 020417 hkh Send optPar back in Result.SOL
% 040103 hkh Revision for v4.2, call iniSolve and endSolve
% 040602 med Help fix PriLev to PriLevOpt
% 041202 hkh If internal differentiation, different call to iniSolve
% 041202 hkh Revise calls to defblbu and StateDef
% 041202 hkh Avoid setting optPar(50) Hessian, unless warm start
% 041222 med Safeguard added for x_0
% 050306 med Step limit added as control parameter
% 050614 med Updated help for optPar vector
% 050614 med Updated iState help, removed optPar 3 and 4
% 050726 med cErrCompute added to avoid incorrect return codes
% 050829 med Removed commented code
% 051216 med optPar(13) default to -1 (no derivative check)
% 060206 med estConsPattern only called if ConsDiff > 6
% 070316 med cErrCompute corrected
% 070613 ang iniSolve call corrected (ConsDiff~=6)
% 080227 hkh Avoid call to SOLSet
% 080604 hkh Safeguard call to iniSolve for empty NumDiff and ConsDiff
% 090827 hkh Set optPar(21)=1E-10 if any pair of finite bounds, n < 1000

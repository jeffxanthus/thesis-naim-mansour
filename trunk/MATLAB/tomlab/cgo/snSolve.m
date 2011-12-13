% snSolve solves the s_n (RBF interpolation surface) subproblem in the RBF 
% based algorithms, of the form:
%
%    min   f(x) = s_n 
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Any of the x could be set as integer valued
%
% Calling syntax:
%
% function snResult = snSolve(snProb)
%
% INPUT PARAMETERS
%
% snProb      Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   sn_f, the RBF interpolation using n sample points in X
%   FUNCS.c   rbf_c, scales x and calls the user constraint function
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   PriLev    Print level in snSolve, >0 gives a printing
%             PriLev-4 used as print level in the 1st (& 2nd) call to tomRun
%             = 1: snSolve cycle and multiMin ExitText output
%             = 1: Error output (if glcCluster), and recovery tomRun output
%             = 2: Number of local steps (if LOCAL)
%             PriLev-4 is used as print level in tomRun local solve call
%             = 4: New best xOpt found during local search (if LOCAL)
%             = 3: Best f(x) from global and local search (if LOCAL)
%             = 3: Best f(x) in global search, globalSolver (if ~LOCAL)
%             = 4: New x found (print of xLoc)
%
%             = 1: onB for all solutions in xOpt
%             = 1: if snSolve picks INTERIOR local minimum instead of global
%                  minimum, print index of local optimum and fOpt(index),
%                  and global optimum denoted: xOpt(:,1) and fOpt(1)
%             = 2: Number of and index of interior points (if more than 1), as
%                  well as f(x) values, called fOpt
%             = 3: Also x values for interior points (if more than 1)
%  
%   c         User constraint function Prob.FUNCS.c
%   dc        User constraint gradient Prob.FUNCS.dc
%   d2c       User constraint Hessian Prob.FUNCS.d2c
%
%   snP       sub point number = Iteration number
%   epsX      Tolerance if vectors in X  are close
%  
% if snProb.LOCAL is true, use
%   MaxIter   Maximal number of iterations used in the local optimization 
%   dLin      Number of linear constraints
% --------------------------------------------
% optParam    Structure in Prob, Prob.optParam 
% ---------------------------------------------
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print one information line each iteration, and the new x tried
%             Default IterPrint = 1. See OUTPUT PRINTING below
%
% ------------------
% Fields in Prob.CGO
% ------------------
%
% globalSolver Solver used for global optimization on the RBF surface
%              If the globalSolver is glcCluster, the fields
%              Prob.GO.maxFunc1, Prob.GO.maxFunc2 and Prob.GO.maxFunc3 are used
%              See the help for maxFunc1, maxFunc2, maxFunc3 in glcCluster
% localSolver  Solver used for local optimization on the RBF surface
%
% X            The n sampled points
%
% SCALE        0-Original search space
%              1-Transform search space to unit cube (default).
% ---------------------------------------------------------
% Fields in Prob.GO (Default values are set for all fields)
% ---------------------------------------------------------
%
% MaxFunc      Maximal number of function evaluations in each global search
% MaxIter      Maximal number of iterations in each global search
% maxFunc1     glcCluster parameter, max function evaluations 1st call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc2     glcCluster parameter, max function evaluations 2nd call
%              Only used if globalSolver is glcCluster, see help globalSolver
% maxFunc3     glcCluster parameter, max sum of function evaluations in 
%              repeated 1st calls trying to get feasible
%              Only used if globalSolver is glcCluster, see help globalSolver
% localSolver  The local solver used by glcCluster
% DIRECT       DIRECT method used in glcCluster: glcSolve or glcDirect(default)
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:  
%             If empty, all variables are assumed non-integer 
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% ---------------------------------------------
%
% snResult    Structure with results from optimization
%  x_k      Matrix with the best points as columns, f(x_k) == f_k.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag Status number
%  Inform   See the snProb.CGO.globalSolver flags
%
% ---------------------------------------
% If globalSolver multiMin is used, the following output is returned:
% ---------------------------------------
% multiMin.xOpt  The matrix of local minima found
% multiMin.fOpt  The k function values in the local optima xOpt(:,i),i=1,...,k
% multiMin.IX    Indices for interior local minima
%
% ---------------------------------------
% For other solvers, the multiMin fields are also defined as:
% ---------------------------------------
% multiMin.xOpt  The global optimum (or global optima)
% multiMin.fOpt  The function value(s) for the global optima
% multiMin.IX    Indices for interior global minima

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 25, 2006. Last modified Oct 20, 2009.

function snResult = snSolve(snProb)

if nargin < 1
   error('snSolve needs one parameter, the structure snProb');
end

% Pick up local parameters
globalSolver          = snProb.CGO.globalSolver;
IntVars               = snProb.MIP.IntVars;
PriLev                = snProb.PriLev;
LOCAL                 = snProb.LOCAL;

X                     = snProb.CGO.X;
F                     = snProb.CGO.F;
nX                    = length(F);
if nX > 100
   % Select only part of sample points in X - now 50 points
   [F0 ix ]           = sort(F);
   X0                 = X(:,ix(1:50));
   F0                 = F0(1:50);
else
   Fmed               = median(F);
   ix                 = find(F <= Fmed);
   X0                 = X(:,ix);
   F0                 = F(ix);
end
k                     = size(X0,2);
epsX                  = snProb.epsX;
if strcmpi(globalSolver,'multiMin') | strcmpi(globalSolver,'multiMINLP')
   snProb.xInit       = X0;
   snProb.fEqTol      = 1E-5;
   snProb.xEqTol      = 1E-4;
   snProb.PriLevOpt   = 0;
elseif strcmpi(globalSolver,'glcCluster')
   snProb.x_0         = [];
   snProb.X0          = X0;
elseif strcmpi(globalSolver,'minlpBB')
   %HKH fix
   if length(IntVars) == snProb.N
      if nX <= 100
         [F000,ix]            = sort(F0);
      end
      % Add some additional initial points beside x_0 = x_min
      snProb.X0               = X0(:,ix(2:min(4,length(ix))));
      snProb.CGO.localSolver  = 'minlpBB';
      snProb.dLin             = snProb.mLin;
      snProb.MaxIter          = 200;
      snProb.optParam.MaxIter = 500;
      LOCAL                   = 1;
   end
end

% HKH? is init of globals necessary - check!
global NLP_x NLP_f NARG
NLP_x=[]; NLP_f=[]; NARG = [];
snResult = tomRunFast(globalSolver,snProb);
PrintResult(snResult,PriLev-4);
modN     = snProb.CGO.modN;


if strcmpi(globalSolver,'multiMin')
   snProb.WarmStart = 1;
   snProb           = WarmDefGLOBAL('multiMin',snProb,snResult);
   snProb.xInit     = max(50,300-k);
   snProb.Mx0       = 20;
   snProb.x_0       = []; % x_0 is part of X, already used in local search
   snResult         = tomRunFast(globalSolver,snProb);
   PrintResult(snResult,PriLev-4);
   snProb.WarmStart = 0;
   if PriLev > 0
      fprintf('snSolve cycle %d: %s ',modN,snResult.ExitText)
      fprintf('x Tolerance %e',snProb.xEqTol);
      fprintf('\n')
   end
end

ExitFlag = snResult.ExitFlag;

if ExitFlag > 0
   if strcmpi(globalSolver,'glcCluster') 
      snProb.WarmStart = 1;
      if snProb.optParam.IterPrint | PriLev > 0
         fprintf('Local min_sn with global solver: No feasible point, ');
         fprintf('Warm Start %s',globalSolver);
         fprintf(' with MaxFunc %d',snProb.optParam.MaxFunc);
         fprintf(' and MaxIter %d\n',snProb.optParam.MaxIter);
      end
      snResult         = tomRunFast(globalSolver,snProb);
      PrintResult(snResult,double(PriLev > 0))
      snProb.WarmStart = 0;
   elseif 0
      % Try the global solver
      mFunc   = snProb.optParam.MaxFunc;
      mIter   = snProb.optParam.MaxIter;
      snProb.optParam.MaxFunc = GOMaxFunc;
      snProb.optParam.MaxFunc = GOMaxIter;
      snResult    = tomRunFast(globalSolver,snProb);
      PrintResult(snResult,double(PriLev > 0))
      snProb.optParam.MaxFunc = mFunc;
      snProb.optParam.MaxIter = mIter;
   end
end

xGlob  = snResult.x_k;
fGlob  = snResult.f_k;

if LOCAL
   localSolver = snProb.CGO.localSolver;
   MaxIter     = snProb.MaxIter;
   dLin        = snProb.dLin;
   iBest       = 0;
   xBest       = [];
   fLoc        = inf;
   i1          = 1;
   if strcmpi(globalSolver,'minlpBB') & length(IntVars) == snProb.N
   %HKH fix
      % No point in doing a search from xGlob(:,1)
      i1          = size(xGlob,2)+1;
      xGlob       = [xGlob,snProb.X0];
      LocSteps    = size(xGlob,2);
      iBest       = 1;
      xBest       = xGlob(:,1);
      fLoc        = fGlob;
   elseif isempty(xGlob)
      % Be sure to start from best point on surface, snProb.x_0, but also in X
      xGlob       = [snProb.x_0,X];
      LocSteps    = min(101,size(xGlob,2));
   else
      % Be sure to start from best point on surface, snProb.x_0
      xGlob       = [xGlob,snProb.x_0];
      LocSteps    = min(21,size(xGlob,2));
   end
   if LocSteps > 1 & PriLev > 1 & size(xGlob,2) > LocSteps
      fprintf('Do %d local search steps ',LocSteps);
      fprintf('out of %d\n',size(xGlob,2));
      fprintf('\n');
   end
   % Do local search from globally best points
   snProb.snP  = 0; % If printout - problem number 0 will be displayed
   snProb.optParam.MaxFunc  = MaxIter*max(snProb.N,10); 
   snProb.optParam.MaxIter  = MaxIter; 
   if ~isempty(IntVars)
      if strcmpi(localSolver,'minlpBB') | strcmpi(localSolver,'oqnlp')
        FixInts = false;
      else
        FixInts = true;
      end
   else
      FixInts = false;
   end
   if FixInts
      x_LL = snProb.x_L(IntVars);
      x_UU = snProb.x_U(IntVars);
   end

   for i = i1:LocSteps
       snProb.x_0  = xGlob(:,i);
       if FixInts
          % Fix integer variables values to nearest integer point
          x00 = max(x_LL,min(x_UU,round(snProb.x_0(IntVars))));
          snProb.x_0(IntVars) = x00;
          snProb.x_L(IntVars) = x00;
          snProb.x_U(IntVars) = x00;
       end
       snR = tomRunFast(localSolver,snProb);
       if PriLev > 4
	  snR.Prob.Name = [snR.Prob.Name,' Local #',num2str(i)];
          PrintResult(snR,PriLev-4)
       elseif PriLev == 4
          fprintf('snSolve: Local  f(x) %30.20f (#%2d of %2d) with %s. ',...
                   snR.f_k,i,LocSteps,localSolver);
          fprintf('CPU %f sec. Iter %d\n', snR.CPUtime,snR.Iter);
       end
       %if snR.CPUtime > 240
       %   snProb.DUNDEE.optPar(1) = 1111111;
       %   snR = tomRunFast(localSolver,snProb);
       %   PrintResult(snR,PriLev-4)
       %   snProb.DUNDEE.optPar(1) = 0;
       %   keyboard
       %end
       if FixInts
          % Reset bounds
          snProb.x_L(IntVars) = x_LL;
          snProb.x_U(IntVars) = x_UU;
       end
       OK  = 1;
       if dLin > 0
          % Check linear constraints
          Ax = snProb.A*snR.x_k;
          AxMax = max(abs(Ax));
          if any(snProb.b_L - Ax > max(1E-5,1E-5*AxMax)  ...
               | Ax - snProb.b_U > max(1E-5,1E-5*AxMax))
             %'DO NOT ACCEPT LOCAL POINT'
             %'LINEAR CONSTRAINTS NOT FULFILLED'
             OK = 0;
          end
       end
       if snR.f_k < fLoc & OK
          iBest    = i;
          snResult = snR;
          xBest    = snR.x_k;
          fLoc     = snR.f_k;
          %if PriLev == 4
          %   xprint(xBest,'xOpt: ');
          %end
       end
   end
   if iBest == 0
      % Local search did not find a better point
      xLoc                     = xGlob(:,1);
      fLoc                     = fGlob;
   else
      xLoc                     = xBest;
   end
   %if LocSteps > 1 & iBest > 1 & PriLev > -2
   %   fprintf('Best local point found in search step %d\n',iBest);
   %end

   %if PriLev > -2 
   if PriLev > 2 
      fprintf('snSolve: Global f(x) %30.20f with %s\n',fGlob,globalSolver);
      fprintf('snSolve: Local  f(x) %30.20f (#%2d out of %3d) with %s\n',...
               fLoc,iBest,LocSteps,localSolver);
   end
else
   if isempty(xGlob)
      xLoc       = [];
      fLoc       = Inf;
   else
      if ~isempty(IntVars)
         x_k  = snResult.x_k(:,1);
         ix = find(abs(x_k(IntVars)-round(x_k(IntVars))) > 1E-10);
         if ~isempty(ix)
            % Fix integer variables values to nearest integer point and rerun
            x00 = max(snProb.x_L(IntVars),min(snProb.x_U(IntVars),...
                      round(x_k(IntVars))));
            snProb.x_0 = x_k;
            snProb.x_0(IntVars(ix)) = x00(ix);
            snProb.x_L(IntVars(ix)) = x00(ix);
            snProb.x_U(IntVars(ix)) = x00(ix);
            snResult = tomRunFast(globalSolver,snProb);
            PrintResult(snResult,PriLev-4)
            % ExitFlag = snResult.ExitFlag;
            % Reset bounds
            % snProb.x_L = x_LL;
            % snProb.x_U = x_UU;
            x_k  = snResult.x_k(:,1);
            fGlob  = snResult.f_k;
         end
         % Make sure we have perfect integer values
         snResult.x_k(IntVars,1) = round(x_k(IntVars)); 
         xGlob  = snResult.x_k;
      end
      xLoc       = xGlob(:,1);
      fLoc       = fGlob;
   end
   if PriLev > 2
      fprintf('Global f(x) %30.20f with %s\n',fGlob,globalSolver);
   end
end
if PriLev > 3
   xprint(xLoc,'x: ');
end

snResult.xLoc = xLoc;
snResult.fLoc = fLoc;

if isempty(snResult.x_k)
   % Should not occur
   if snResult.ExitFlag == 0
      fprintf('snSolve - local problem: Failure in global solver!')
      fprintf(' Not feasible\n')
   else
      fprintf('snSolve - local problem: Failure in global solver!\n')
   end
   fprintf('ExitFlag = %d. ExitText: ', snResult.ExitFlag);
   fprintf('%s.', snResult.ExitText);   
   snResult.f_k = Inf;
end

if strcmpi(globalSolver,'multiMin')
   nLocal                 = snResult.multiMin.Info.nLocal;
   nGlobal                = snResult.multiMin.Info.nGlobal;
   xOpt                   = snResult.multiMin.XX(:,1:nLocal);
   fOpt                   = snResult.multiMin.FX(1,1:nLocal);
   Inform                 = snResult.multiMin.EX(3,1:nLocal);
   k                      = 1;
else
   nLocal                 = 1;
   nGlobal                = 1;
   xOpt                   = snResult.x_k;
   fOpt                   = snResult.f_k;
   k                      = size(xOpt,2);
   %snResult.multiMin.xOpt = xOpt;
   %snResult.multiMin.fOpt = fOpt;
end
snResult.multiMin.nLocal  = nLocal;
snResult.multiMin.nGlobal = nGlobal;
snResult.multiMin.xOpt    = xOpt;
snResult.multiMin.fOpt    = fOpt;
onB                       = nOnBound(xOpt,snProb.x_L,snProb.x_U,epsX);
IX                        = find(onB == 0);
snResult.multiMin.onB     = onB;
snResult.multiMin.IX      = IX;
%HKH fix
if PriLev > 1
   xprinti(onB,'   onB:');
end
if size(xOpt,2) == 1
   snResult.x_k   = xOpt;
   snResult.f_k   = fOpt(1);
else
   if nGlobal == 1 
      snResult.x_k   = xOpt;
      snResult.f_k   = fOpt(1);
   else
      oB = min(onB(1:nGlobal));
      iB      = find(onB == oB);
%HKH fix
      if PriLev > 1
         fprintf('snSolve: oB %d ',oB);
         fprintf('iB %d ',iB);
      end
      if length(iB) == 1
         snResult.x_k   = xOpt(:,iB);
         snResult.f_k   = fOpt(iB);
      else
         % Select global optimum with lowest components on bound
         % and at largest distance from sampled points X
         [Dist iD]      = max(min(tomsol(30,xOpt(:,iB),X)'));
         ix = iB(iD);
%HKH fix
         if PriLev > 1
            fprintf('   snSolve: Dist %f ',Dist);
            fprintf('iD %d ',iD);
            fprintf('ix %d ',ix);
            fprintf('fOpt-1 %16.10f ',fOpt(1));
            fprintf('fOpt-ix %16.10f ',fOpt(ix));
         end
         snResult.x_k   = xOpt(:,ix);
         snResult.f_k   = fOpt(ix);
      end
%HKH fix
      if PriLev > 1
         fprintf('\n')
      end
   end
end
if 0
if ~isempty(IX)
   if k > 1 & IX(1) ~= 1
      snResult.x_k   = xOpt(:,IX(1)); % Select interior point, if possible
      snResult.f_k   = fOpt(min(length(fOpt),IX(1)));
      if PriLev > 0 
         fprintf('snSolve: *** PICK Interior local minimum on RBF surface: ')
         fprintf(' %d',IX)
         fprintf(' fOpt(1) %f',fOpt(1))
         fprintf(' fOpt(%d) %f',IX(1),snResult.f_k)
         fprintf('\n')
	 xprint(xOpt(:,1),'xOpt(:,1):');
      end
   end
elseif size(xOpt,2) > 1
   % Select only one of the x solutions, multiMINLP may give plenty
   % Could select the one furthest away from X
   snResult.x_k   = xOpt(:,1);
   snResult.f_k   = fOpt(1);
end
end

if PriLev > 0 
   xprinti(onB(1:min(50,length(onB))),'onB:');
   if PriLev > 1 & length(IX) > 1
      fprintf('snSolve: %d interior points !!! IX:',length(IX))
      fprintf(' %d',IX)
      if length(IX) > 10, fprintf('\n'); end
      fprintf(' fOpt:')
      fprintf(' %f',fOpt(IX))
      fprintf('\n')
      if PriLev > 2 & modN == snProb.CGO.N
         disp(xOpt(:,IX))
      end
   end
end
% fprintf('snSolve  f_k = %20.10f\n',snResult.f_k);

% MODIFICATION LOG
%
% 061125  hkh  First version written from code in rbfSolve, revised, simplified
% 070223  hkh  Use snProb.x_L and snProb.x_U, not just x_L, x_U
% 070224  hkh  IterPrint must be picked from snProb.optParam.IterPrint
% 070302  hkh  Revise comments, clean up
% 071002  hkh  LOCAL wrong for multiMin. Use snProb.LOCAL instead.
% 071005  hkh  Use both X and random initial points running multiMin
% 071009  hkh  Use snProb.x_0 in LOCAL search if nonempty
% 071010  hkh  Use snProb.CGO.X in glcCluster, input snProb.X0
% 080410  hkh  Wrong reset in snProb.x_L/x_U if IntVars
% 080414  hkh  Incorrect code in local search
% 080608  hkh  Use tomRunFast
% 080623  hkh  Call PrintResult separately
% 080626  hkh  Limit onB output to 50 elements
% 080626  hkh  Utilize minlpBB with local searches for pure IP problems
% 090824  hkh  Define global NLP_x NLP_f NARG as empty before tomRunFast
% 090911  hkh  Only return one xOpt vector. Add test for multiMINLP
% 090911  hkh  Avoid crash accessing fOpt(IX(1)), if many min from multiMINLP
% 090918  hkh  Select X with Fpen <= median(Fpen) if n <= 100, otherwise
%              if n > 100, use the 50 points from  X with lowest Fpen
% 091020  hkh  New use of revised multiMin. Make better choice of best global x_k

% gnSolve solves the gn (Target value) subproblem in the RBF based algorithms, 
% of the form:
%
%    min   f(x) = g_n
%     x
%    s/t   x_L <= x <= x_U, x_L and x_U finite
%          b_L <= A x  <= b_U
%          c_L <= c(x) <= c_U
%
% Any of the x could be set as integer valued
%
% Calling syntax:
%
% function gnResult = gnSolve(gnProb)
%
% INPUT PARAMETERS
%
% gnProb      Structure, where the following variables are used:
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
%   X0        A matrix of initial search points in local search
%   PriLev    Print level in gnSolve, == 0 is silent, except if failure
%             PriLev-4 used as print level in the 1st (& 2nd) call to tomRun
%             = 1: gnSolve modN cycle and multiMin ExitText output
%             = 3: Dmin - distance from every solution to closest in X
%             = 1: Header text and indices for xOpt solutions too close to X
%             = 3: Up to 30 of the deleted xOpt solutions, too close to X
%             = 1: Best pnt 1, f(x) and x value (up to 10 components), Inform
%             = 2: The up to 15 best f(x) values with multiMin
%             = 2: fOpt, modN, xOpt(1:10,:) for <= 5 best solutions
%             = 3: Compute, print my,Sq = (sn_f-fnStar)^2, GN_f in the same row 
%                  if isinf(fnStar), sn_f is printed, GN_f = -1/my (modN == -1)
%                  if ~isinf(fnStar), GN_f = -1/(my*(sn_f-fnStar)^2)
%             = 1: onB for the up to 50 best solutions in xOpt
%             = 1: List of indices for the interior points found
%             = 2: Indices, fOpt, xOpt(1:10,:),Inform for <= 5 best interiors
%             = 4: Best f(x) from global search
%             = 2: Number of local steps (if LOCAL)
%             = 5: Print local search number  (if LOCAL)
%             PriLev-4 is used as print level in tomRun local solve call
%             = 4: x_0 in local search (if LOCAL)
%             = 4: New best xOpt found during local search (if LOCAL)
%             = 3: Best f(x) from global and local search (if LOCAL)
%             = 3: Best f(x) in global search, globalSolver (if ~LOCAL)
%             = 4: New x found (print of xLoc)
%  
%   c         User constraint function Prob.FUNCS.c
%   dc        User constraint gradient Prob.FUNCS.dc
%   d2c       User constraint Hessian Prob.FUNCS.d2c
%
%   LOCAL
%   IX        (only if multiMin, iteration dependent, output from snSolve)
%   xOpt      (only if multiMin, iteration dependent, output from snSolve)
%   fOpt      (only if multiMin, iteration dependent, output from snSolve)
%   epsX      (only if multiMin)
%  
% if gnProb.LOCAL is true, use
%   MaxIter   Maximal number of iterations used in the local optimization 
%   dLin      Number of linear constraints
%
% ------------------
% Fields in Prob.CGO
% ------------------
%
% modN         Cycle step
%
% fnStar       0 (=inf) or target value fnStar
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
%
% TargetMin    Which minimum of several to pick in target value problem
%              =0 Use global minimum
%              =1 Use best interior local minima, if none use global minimum
%              =2 Use best minimum with lowest number of coefficients on bounds
%
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
% OUTPUT PARAMETERS
% ---------------------------------------------
%
% gnResult    Structure with results from optimization
%  x_k      Matrix with the best points as columns, f(x_k) == f_k.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag Status number
%  Inform   See the gnProb.CGO.globalSolver flags
%
% ---------------------------------------
% If globalSolver multiMin is used, the following output is returned:
% ---------------------------------------
% multiMin.xOpt  The matrix of local minima found
% multiMin.fOpt  The k function values in the local optima xOpt(:,i),i=1,...,k
% multiMin.IX    Indices for interior local minima

%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 24, 2006. Last modified Oct 21, 2009.

function gnResult = gnSolve(gnProb)

if nargin < 1
   error('gnSolve needs one parameter, the structure gnProb');
end
global NLP_x NLP_f NARG % Used if PriLev > 2

% Pick up local parameters
%globalSolver = 'oqnlp'
globalSolver  = gnProb.CGO.globalSolver;
IntVars       = gnProb.MIP.IntVars;
PriLev        = gnProb.PriLev;
LOCAL         = gnProb.LOCAL;
TargetMin     = gnProb.CGO.TargetMin;

if strcmpi(globalSolver,'multiMin') | strcmpi(globalSolver,'multiMINLP')
   epsX             = gnProb.epsX;
   gnProb.xInit     = gnProb.X0;
   gnProb.xEqTol    = 1E-4;
   gnProb.fEqTol    = 1E-5;
   gnProb.PriLevOpt = 0;
elseif strcmpi(globalSolver,'glcCluster')
   % In gnProb.X0 is already the matrix of extra initial points to be used
elseif strcmpi(globalSolver,'minlpBB')
   %HKH fix
   if length(IntVars) == gnProb.N
      gnProb.MIP.IntVars      = [];
      gnProb.xInit            = 100;
      gnProb.xEqTol           = 1E-1;
      PriLevOpt               = gnProb.PriLevOpt;
      gnProb.PriLevOpt        = 0;
      gnProb.PriLevOpt        = 1;
      gnResult                = tomRunFast('multiMin',gnProb);
      PrintResult(gnResult, PriLev-4);
      %xprint(gnResult.multiMin.fOpt);
      k                       = min(5,size(gnResult.multiMin.xOpt,2));
      gnProb.X0               = gnResult.multiMin.xOpt(:,1:k);
      gnProb.PriLevOpt        = PriLevOpt;
      gnProb.MIP.IntVars      = IntVars;
      gnProb.CGO.localSolver  = 'minlpBB';
      gnProb.dLin             = gnProb.mLin;
      gnProb.MaxIter          = 200;
      gnProb.optParam.MaxIter = 500;
      LOCAL                   = 1;
   end
end
gnResult = tomRunFast(globalSolver,gnProb);
PrintResult(gnResult, PriLev-4);
%if gnResult.CPUtime > 240
%   gnProb.DUNDEE.optPar(1) = 1111111;
%   gnResult = tomRunFast(globalSolver,gnProb);
%   PrintResult(gnResult,PriLev-4)
%   gnProb.DUNDEE.optPar(1) = 0;
%   'global'
%   keyboard
%end
if strcmpi(globalSolver,'multiMin')
   gnProb.WarmStart          = 1;
   %gnProb.xInit              = max(50,300-size(gnProb.X0,2));
   gnProb.xInit              = max(50,100-size(gnProb.X0,2));
   %gnProb.Mx0                = ceil(0.95*gnProb.xInit);
   gnProb.Mx0                = ceil(0.90*gnProb.xInit);
   gnProb                    = WarmDefGLOBAL('multiMin',gnProb,gnResult);
   gnResult                  = tomRunFast(globalSolver,gnProb);
   PrintResult(gnResult, PriLev-4);
   gnProb.WarmStart          = 0;
   nLocal                    = gnResult.multiMin.Info.nLocal;
   nGlobal                   = gnResult.multiMin.Info.nGlobal;
   nSolution                 = gnResult.multiMin.Info.nSolution;
   if nLocal > 0
      xOpt                     = gnResult.multiMin.XX(:,1:nLocal);
      fOpt                     = gnResult.multiMin.FX(1,1:nLocal);
      Inform                   = gnResult.multiMin.EX(3,1:nLocal);
   else
      xOpt                     = gnResult.multiMin.XX(:,1:nSolution);
      fOpt                     = gnResult.multiMin.FX(1,1:nSolution);
      Inform                   = gnResult.multiMin.EX(3,1:nSolution);
   end
   gnResult.multiMin.nLocal    = nLocal;
   gnResult.multiMin.nGlobal   = nGlobal;
   gnResult.multiMin.nSolution = nSolution;
   gnResult.multiMin.xOpt      = xOpt;
   gnResult.multiMin.fOpt      = fOpt;
   gnResult.multiMin.Inform    = Inform;
   modN                        = gnProb.CGO.modN;
   if PriLev > 0
      fprintf('gnSolve cycle %d: %s ',modN,gnResult.ExitText)
      fprintf('x Tolerance %e',gnProb.xEqTol);
      fprintf('\n')
   end
   Dist   = tomsol(30,xOpt,gnProb.CGO.X);
   Dmin   = min(Dist');
%fix
   if PriLev > 2
      xprint(Dmin(1:min(length(Dmin),15)),'   Dmin:',[],15);
   end
   ix     = find(Dmin > 1E-6);
   if length(ix) ~= size(xOpt,2);
      ixD = find(Dmin <= 1E-6);
      if PriLev > 0
         fprintf('****** Delete solutions in xOpt, too close to X\n');
         xprinti(ixD(1:min(length(ixD),30)),'idx: ');
         if PriLev > 2
            for ii=1:min(30,length(ixD))
	        xprint(xOpt(:,ixD(ii)),'xOpt: ');
            end
         end
      end
      gnResult.multiMin.xOpt   = gnResult.multiMin.xOpt(:,ix); 
      gnResult.multiMin.fOpt   = gnResult.multiMin.fOpt(ix); 
      gnResult.multiMin.Inform = gnResult.multiMin.Inform(ix); 
      xOpt                     = gnResult.multiMin.xOpt;
      fOpt                     = gnResult.multiMin.fOpt;
      Inform                   = gnResult.multiMin.Inform;
   end
   if PriLev > 0
      fprintf('gnSolve: Best pnt # %2d: ',1)
      fprintf('f(x) %18.14f ',fOpt(1))
      fprintf('x: ')
      fprintf('%9.7f ',xOpt(1:min(10,end),1))
      fprintf(' Inform %d',Inform(1))
      fprintf('\n')
      if PriLev > 1
         xprint(fOpt,'fOpt:',[],15);
      end
   end
   k                     = size(xOpt,2);
   onB                   = nOnBound(xOpt,gnProb.x_L,gnProb.x_U,epsX);
   IX                    = find(onB == 0);
   gnResult.multiMin.onB = onB;
   gnResult.multiMin.IX  = IX;
   if PriLev > 1
      fnStar = gnProb.CGO.fnStar;
      for i = 1:min(5,length(fOpt))
          fprintf('# %2d: ',i)
          fprintf('f(x) %14.10f ',fOpt(i))
          if PriLev > 2
             gnProb.CGO.modN = -1;
             NLP_x=[]; NLP_f=[]; NARG = [];
             my = -1/gn_f(xOpt(:,i),gnProb);
             gnProb.CGO.modN = modN;
             fprintf('my %14.10f ',my)
             % Compute sn_f
	     if modN == -1
                snf = tomsol(21, xOpt(:,i));
                fprintf('snf %14.10f ',snf)
             else 
                Sq = (tomsol(21, xOpt(:,i))-fnStar)^2;
                fprintf('Sq  %14.10f ',Sq)
             end 
             %fprintf('f %14.10f ',-log(1/(f1*f2)))
	     %if modN == -1
             %   fprintf('GN_f %14.10f ',-1/my)
             %   fprintf('f %9.3f ',fOpt(i)/(-1/my))
             %else 
             %   fprintf('GN_f %14.10f ',-1/(Sq*my))
             %   fprintf('f %9.3f ',fOpt(i)/(-1/(Sq*my)))
             %end 
          end 
          fprintf('modN %d ',modN)
          fprintf('x: ')
          fprintf('%9.7f ',xOpt(1:min(10,end),i))
          fprintf('\n')
      end
   end
%fix
   if PriLev > 0
      xprinti(onB,'   onB:');
   end
   if PriLev > 0 & ~isempty(IX)
      fprintf('gnSolve: Interior points: ')
      fprintf(' %d',IX)
      fprintf('\n')
      if PriLev > 1
         for i = 1:min(5,length(IX))
             fprintf('# %2d: ',IX(i))
             fprintf('f(x) %18.14f ',fOpt(IX(i)))
             fprintf('x: ')
             fprintf('%9.7f ',xOpt(1:min(10,end),IX(i)))
             fprintf(' Inform %d',Inform(IX(i)))
             fprintf('\n')
         end
      end
      %if IX(1) > 1
      %   fprintf('Use interior GN surface point\n')
      %   gnResult.x_k = xOpt(:,IX(1));
      %end
   end     
   if (TargetMin == 1 | TargetMin == 2) & (isempty(IX) | IX(1) ~= 1) 
      if ~isempty(IX)
         % Use best interior point from GN problem
         if IX(1) > 1
            if PriLev > 0
               fprintf('Use interior GN surface point\n')
            end
            gnResult.x_k = xOpt(:,IX(1));
            gnResult.f_k = fOpt(IX(1));
         end
      elseif TargetMin == 2 & modN ~= -1
         IXsn   = gnProb.IX;
         if isempty(IX) & length(IXsn) == 2
            % Use second best interior point from RBF surface
            xOptsn = gnProb.xOpt;
            fOptsn = gnProb.fOpt;
            if PriLev > 0
               fprintf('Use interior RBF surface point\n')
            end
            gnResult.x_k = xOptsn(:,IXsn(2));
            gnResult.f_k = fOptsn(IXsn(2));
         elseif isempty(IX) & length(IXsn) > 2
            % Use interior point from RBF surface most far from X
            xOptsn = gnProb.xOpt;
            fOptsn = gnProb.fOpt;
            xxD = zeros(length(IXsn),1);
            X = gnProb.CGO.X;
            for ii=2:length(IXsn)
                xxD(ii) = min(tomsol(30,xOptsn(:,IXsn(ii)),X));
            end
            [mm ii]=max(xxD);
            gnResult.x_k = xOptsn(:,IXsn(ii));
            gnResult.f_k = fOptsn(IXsn(ii));
            if PriLev > 0
               fprintf('Use interior RBF surface point WITH MAX DIST, ')
               fprintf('Pnt %d\n',ii)
            end
         end
      end
   elseif TargetMin == 3
      STRICT = 0;
      if STRICT == 1
         [onBMin, onIX] = min(onB(1:nGlobal));
      else
         [onBMin, onIX] = min(onB);
      end
      if PriLev > 0
         %fix
         fprintf('   gnSolve: onBMin %d ',onBMin);
         fprintf('onIX: ');
         fprintf('%d ',onIX);
         if length(fOpt) >= onIX
            fprintf('fOpt-1 %16.10f ',fOpt(1));
            fprintf('fOpt-ix %16.10f ',fOpt(onIX));
         else
            fprintf('length(fOpt) %d\n',length(fOpt));
         end
         fprintf('\n')
      end
      if PriLev > 0 & onIX ~= 1
         fprintf('Use GN minimum %d with lowest onB = %d\n',onIX,onBMin)
      end
      gnResult.x_k = xOpt(:,onIX);
      gnResult.f_k = fOpt(onIX);
   end
elseif strcmpi(globalSolver,'glcCluster')
   gnProb.WarmStart = 1;
   for i = 1:5
       if isempty(gnResult.x_k)
          fprintf('gnSolve: Warm Start %d. glcCluster was not feasible\n',i)
          gnResult = tomRunFast(globalSolver,gnProb);
       end
       PrintResult(gnResult, PriLev-4);
   end
   gnProb.WarmStart = 0;
elseif strcmpi(globalSolver,'glcDirect') | strcmpi(globalSolver,'glcFast')
   gnProb.WarmStart = 1;
   for i = 1:5
       if isempty(gnResult.x_k)
          fprintf('gnSolve: Warm Start %d. ',i)
          fprintf('%s was not feasible\n',globalSolver)
          gnResult = tomRunFast(globalSolver,gnProb);
       end
       PrintResult(gnResult, PriLev-4);
   end
   gnProb.WarmStart = 0;
elseif strcmpi(globalSolver,'glcSolve') | strcmpi(globalSolver,'glbSolve')
   gnProb.WarmStart = 1;
   for i = 1:2
       if isempty(gnResult.x_k)
          fprintf('gnSolve: Warm Start %d. ',i)
          fprintf('%s was not feasible\n',globalSolver)
          gnResult = tomRunFast(globalSolver,gnProb);
       end
       PrintResult(gnResult, PriLev-4);
   end
   gnProb.WarmStart = 0;
end

if isempty(gnResult.x_k)
   % Should never occur
   if gnResult.ExitFlag == 0
      fprintf('gnSolve - global problem: Failure in global solver!')
      fprintf(' Not feasible\n')
   else
      fprintf('gnSolve - local problem: Failure in global solver!\n')
   end
   fprintf('ExitFlag = %d. ExitText: ', gnResult.ExitFlag);
   fprintf('%s.', gnResult.ExitText);   
%HKH fix now when multiMin may return empty result
   %backupSolver1 = gnProb.backupSolver1;
   backupSolver1 = 'glcCluster';
   if ~isempty(backupSolver1)
      fprintf('Rerun with %s\n',backupSolver1);
      if strcmpi(backupSolver1,'oqnlp')
         gnProb.OQNLP.options.RANDOMSEED = 12345;
      end
      gnResult = tomRunFast(globalSolver,gnProb);
      PrintResult(gnResult, PriLev-4);
   else
      fprintf('\n')
   end
   if isempty(gnResult.x_k)
      fprintf('WARNING: Still not global feasible. ');
      backupSolver2 = gnProb.backupSolver2;
      if ~isempty(backupSolver2)
         fprintf('Rerun with %s\n',backupSolver2);
         if strcmpi(backupSolver2,'oqnlp')
            gnProb.OQNLP.options.RANDOMSEED = 54321;
         end
         gnResult = tomRunFast(globalSolver,gnProb);
         PrintResult(gnResult, PriLev-4);
      else
	 fprintf('\n')
      end
   end
   if isempty(gnResult.x_k)
      fprintf('WARNING: Complete failure in gnSolve.\n');
      LOCAL = 0; % No point in running the LOCAL block - complete failure
   end
end
xGlob  = gnResult.x_k;
fGlob  = gnResult.f_k;
if PriLev > 3
   if ~isempty(xGlob), xprint(xGlob(:,1),'x: '); end
end
if LOCAL
   localSolver = gnProb.CGO.localSolver;
   MaxIter     = gnProb.MaxIter;
   dLin        = gnProb.dLin;
   maxLocTry   = max(21,gnProb.N*3);
   iBest       = 0;
   xBest       = [];
   fLoc        = inf;
   i1          = 1;
   if strcmpi(localSolver,'minlpBB')
      % No point in doing a search from xGlob(:,1)
      i1          = size(xGlob,2)+1;
      xGlob       = [xGlob,gnProb.X0];
      iBest       = 1;
      xBest       = xGlob(:,1);
      fLoc        = fGlob;
   elseif isempty(gnProb.x_0)
      if ~isempty(gnProb.X0)
         xGlob       = [xGlob,gnProb.X0];
      end
   elseif size(xGlob,2) > maxLocTry-2
      xGlob       = [xGlob(:,1:maxLocTry-1),gnProb.x_0];
   else
      xGlob       = [xGlob,gnProb.x_0];
      if ~isempty(gnProb.X0)
         xGlob       = [xGlob,gnProb.X0];
      end
   end
   LocSteps    = min(maxLocTry,size(xGlob,2));
   %if LocSteps > 1 & PriLev > 1 
   if LocSteps > 1 & PriLev > 1 & size(xGlob,2) > LocSteps
      fprintf('Do %d local search steps ',LocSteps);
      fprintf('out of %d\n',size(xGlob,2));
      fprintf('\n');
   end
   % Do local search from globally best points
   gnProb.snP  = 0; % If printout - problem number 0 will be displayed
   gnProb.optParam.MaxFunc  = MaxIter*max(gnProb.N,10); 
   gnProb.optParam.MaxIter  = MaxIter; 
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
      x_LL = gnProb.x_L(IntVars);
      x_UU = gnProb.x_U(IntVars);
   end

   for i = i1:LocSteps
       gnProb.x_0  = xGlob(:,i);
       if FixInts
          % Fix integer variables values to nearest integer point
          x00 = max(x_LL,min(x_UU,round(gnProb.x_0(IntVars))));
          gnProb.x_0(IntVars) = x00;
          gnProb.x_L(IntVars) = x00;
          gnProb.x_U(IntVars) = x00;
       end
       if PriLev == 7
          xprint(gnProb.x_0,'x0:   ');
       end
       gnR = tomRunFast(localSolver,gnProb);
       if PriLev > 4
	  gnR.Prob.Name = [gnR.Prob.Name,' Local #',num2str(i)];
          PrintResult(gnR, PriLev-4);
       elseif PriLev == 4
          fprintf('gnSolve: Local  f(x) %30.20f (#%2d of %3d) with %s. ',...
                   gnR.f_k,i,LocSteps,localSolver);
          fprintf('CPU %f sec. Iter %d\n', gnR.CPUtime,gnR.Iter);
       end
       %if gnR.CPUtime > 240
       %   gnProb.DUNDEE.optPar(1) = 1111111;
       %   gnR = tomRunFast(localSolver,gnProb);
       %   PrintResult(gnR,PriLev-4)
       %   gnProb.DUNDEE.optPar(1) = 0;
       %   keyboard
       %end
       if FixInts
          % Reset bounds
          gnProb.x_L(IntVars) = x_LL;
          gnProb.x_U(IntVars) = x_UU;
       end
       OK  = 1;
       %if gnR.ExitFlag ~= 0 & dLin > 0
       if dLin > 0
          % Check linear constraints
          Ax = gnProb.A*gnR.x_k;
          AxMax = max(abs(Ax));
          if any(gnProb.b_L - Ax > max(1E-5,1E-5*AxMax)  ...
               | Ax - gnProb.b_U > max(1E-5,1E-5*AxMax))
             %'DO NOT ACCEPT LOCAL POINT'
             %'LINEAR CONSTRAINTS NOT FULFILLED'
             %ExitFlag = gnR.ExitFlag
             %Inform   = gnR.Inform
             OK = 0;
          end
       end
       if gnR.f_k < fLoc & OK
          iBest    = i;
          gnResult = gnR;
          xBest    = gnR.x_k;
          fLoc     = gnR.f_k;
          if PriLev == 7
             xprint(xBest,'xOpt: ');
          end
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

   %if PriLev > -1 
   if PriLev > 2 
      fprintf('gnSolve: Global f(x) %30.20f with %s\n',fGlob,globalSolver);
      fprintf('gnSolve: Local  f(x) %30.20f (#%2d of %3d) with %s\n',...
               fLoc,iBest,LocSteps,localSolver);
   end
   if fLoc > fGlob & gnProb.mNonLin == 0
      xDiff1 = sum(max(0,gnProb.x_L-xGlob(:,1))+max(0,xGlob(:,1)-gnProb.x_U));
      fprintf('gnSolve: WARNING! ');
      fprintf('Global f(x) %30.20f with %s (%e)',fGlob,globalSolver,xDiff1);
      xDiff2 = sum(max(0,gnProb.x_L-xLoc)+max(0,xLoc-gnProb.x_U));
      fprintf(' > Local f(x) %30.20f (#%2d) with %s (%e)\n', ...
              fLoc,iBest,localSolver,xDiff2);
      %if fLoc > fGlob, keyboard, end
      if xDiff1 <= xDiff2
         % Obviously worse solution from local search, use global solution
         xLoc                     = xGlob(:,1);
         fLoc                     = fGlob;
      end
   end
else
   if isempty(xGlob)
      xLoc       = [];
      fLoc       = Inf;
   else
      if ~isempty(IntVars)
         % HKH - really needed with this safe guard?
         x_k  = gnResult.x_k(:,1);
         ix = find(abs(x_k(IntVars)-round(x_k(IntVars))) > 1E-10);
         if ~isempty(ix)
            % Fix integer variables values to nearest integer point and rerun
            x00 = max(gnProb.x_L(IntVars),min(gnProb.x_U(IntVars),...
                      round(x_k(IntVars))));
            gnProb.x_0 = x_k;
            gnProb.x_0(IntVars(ix)) = x00(ix);
            gnProb.x_L(IntVars(ix)) = x00(ix);
            gnProb.x_U(IntVars(ix)) = x00(ix);
            gnResult = tomRunFast(globalSolver,gnProb);
            PrintResult(gnResult, PriLev-4);
            % ExitFlag = gnResult.ExitFlag;
            % Reset bounds
            % gnProb.x_L = x_LL;
            % gnProb.x_U = x_UU;
            x_k    = gnResult.x_k(:,1);
            fGlob  = gnResult.f_k;
         end
         % Make sure we have perfect integer values
         gnResult.x_k(IntVars,1) = round(x_k(IntVars)); 
         xGlob  = gnResult.x_k;
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

gnResult.xLoc = xLoc;
gnResult.fLoc = fLoc;

% fprintf('gnSolve fLoc = %20.10f fGlob = %20.10f\n',fLoc,fGlob);

% MODIFICATION LOG
%
% 061124  hkh  First version written from code in rbfSolve, revised, simplified
% 070302  hkh  Revise comments, clean up
% 071005  hkh  Use both X and random initial points running multiMin
% 071005  hkh  Revised comments and printing
% 071005  hkh  Avoid using worse interior points as final results
% 071008  hkh  New input TargetMin when using multiMin
% 071009  hkh  Use gnProb.x_0 in LOCAL search if nonempty
% 071010  hkh  Use gnProb.X0 in glcCluster and multiMin, and if LOCAL
% 071128  hkh  X=gnProb.CGO.X had to be set at two places
% 080414  hkh  IntVars not defined and used in LOCAL search
% 080414  hkh  IntVars not safe guarded when not LOCAL set
% 080414  hkh  Incorrect code when doing local search
% 080416  hkh  Safe guard, use global solution if local is worse
% 080608  hkh  Use tomRunFast instead of tomRun
% 080623  hkh  Call PrintResult separately 
% 080625  hkh  Utilize minlpBB for pure IP problems
% 090824  hkh  Minor mlint revision
% 090911  hkh  Add test for solver multiMINLP

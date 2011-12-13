% multiMINLP.m
%
% multiMINLP solves general constrained mixed-integer global nonlinear
% optimization problems.
%
% It is aimed for problems where the number of integer combinations nMax
% is huge and relaxation of the integer constraint is possible.
%
% If no integer variables, multiMINLP calls multiMin.
% If nMax <= min(Prob.optParam.MaxFunc,5000), glcDirect is used.
% Otherwise, multiMINLP first finds a set M of local minima calling 
% multiMin with no integer restriction on any variable. The best local 
% minimum is selected. glcDirect is called to find the best integer 
% feasible solution fIP in a small area (< +- 2 units) around the best 
% local minimum found.
%
% The other local minima are pruned, if fOpt(i) > fIP, no integer feasible
% solution could be found close to this local minimum i.
%
% The area close to the remaining candidate local minima are searched one 
% by one by calling glcDirect to find any fIPi < fIP.
%
% multiMINLP solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Usage: See below
%
% Calling syntax:
%
% function Result = tomRun('multiMINLP',Prob)
%
% INPUT PARAMETERS
%
%
% Prob    Structure, defined as to solve a standard MINLP problem. The Prob 
%         structure is fed to the localSolver. See e.g. minlpAssign.
%
%         See multiMin and glcDirect for input to the subsolvers e.g. 
%         Prob.xInit is used in multiMin (and fCut, RandState, xEQTol).
%
%   x_L       Lower bounds for each element in x. If generating random
%             points, -inf elements of x_L are set to min(-L,xMin,x_U-L)
%             xMin is the minimum of the finite x_L values.
%   x_U       Upper bounds for each element in x. If generating random points,
%             inf elements of x_U are set to max(L,xMax,x_L+L)
%             xMax is the maximum of the finite x_U values.
%
%             L is 100 for nonlinear least squares, otherwise 1000.
%
%   b_L       Lower bounds for the linear constraints.
%   b_U       Upper bounds for the linear constraints.
%   A         Linear constraint matrix.
%   c_L       Lower bounds for the nonlinear constraints.
%   c_U       Upper bounds for the nonlinear constraints.
%
%   PriLev    Print Level
%             0 = Silent
%             1 = Display 2 summary rows
%             2 = Display some extra summary rows
%             5 = Print level 1 in tomRun call
%             6 = Print level 2 in tomRun call
%             7 = Print level 3 in tomRun call
%
%   xInit     Used in multiMin. See help for multiMin.
%
% GO
%   localSolver The local solver used to run all local optimizations.
%               Default is the license dependent output of
%               GetSolver('con',1,0).
%
% optParam    Structure in Prob, Prob.optParam. Defines optimization
%             parameters. Fields used:
%   MaxFunc   Max number of function evaluations in each subproblem
%   fGoal     Goal for function value f(x), if empty not used
%             If fGoal is reached, no further local optimizations are done
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal ~= 0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%             Default 1e-8.
%   bTol      Linear constraint feasibility tolerance. Default 1e-6
%   cTol      Nonlinear constraint feasibility tolerance. Default 1e-6
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:
%             If empty, all variables are assumed non-integer
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer
%             variables
%
%   nMax      Number of integer combinations possible, if empty multiMINLP computes nMax
%
%   Rfac      Reduction factor for real variables in MINLP subproblem close to local multiMINLP minimum
%             Bounds set to x_L = max(x_L,x-Rfac*(x_U-x_L)) and x_U = min(x_U,x+Rfac*(x_U-x_L)). Def 0.25
%
% OUTPUT PARAMETERS:
%
% Result    Result structure from the last good optimization step giving the
%           best f(x) value, the possible global MINLP minimum
%
% The following fields in Result are changed by multiMINLP before return:
%
% ExitFlag  = 0 normal output, of if fGoal set and found
%           = 1 A solution reaching the user defined fGoal was not found
%           = 2 Unbounded problem
%           = 4 Infeasible problem
%
% The Solver, SolverAlgorithm and ExitText fields are also reset
%
% A special field in Result is also returned, Result.multiMINLP:
%
%  xOpt     Prob.N x k matrix with k distinct local optima, the test being
%              norm(x_k-xOpt(:,i)) <= xEqTol*max(1,norm(x_k))
%           that if fulfilled assumes x_k to be to close to xOpt(:,i)
%  fOpt     The k function values in the local optima xOpt(:,i),i=1,...,k
%  Inform   The Inform value returned by the local solver when finding
%           each of the local optima xOpt(:,i); i=1,...,k
%           The Inform value can be used to judge the validity of the local
%           minimum reported.
%  localTry Total number of local searches
%  Iter     Total number of iterations
%  FuncEv   Total number of function evaluations
%  GradEv   Total number of gradient evaluations
%  HessEv   Total number of Hessian evaluations
%  ConstrEv Total number of constraint function evaluations
%  ExitText Text string giving Inform information
%
% USAGE:
%
% Driver call, including printing with level 2:
%    Result = tomRun('multiMINLP',Prob,2);
%
% Direct solver call
%    Result = multiMINLP(Prob);
%
% Then print the result of the best local search with the call
%    PrintResult(Result);

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc., $Release: 7.5.0$
% Written Sept 11, 2006.    Last modified May 5, 2010.

function Result = multiMINLP(Prob)

if nargin < 1
   error('multiMINLP needs input structure Prob');
end

probType        = Prob.probType;
CLS             = 0;
if isfield(Prob.GO,'localSolver')
    localSolver = deblank(Prob.GO.localSolver);
else
    localSolver = [];
end
if isempty(localSolver) & probType==checkType('cls')
    SolvType    = checkType('cls');
    localSolver = GetSolver('cls',0,0);
    CLS         = 1;
else
    SolvType    = checkType('con');
end
if isempty(localSolver)
    localSolver = GetSolver('con',1,0);
end
Prob            = ProbCheck(Prob,localSolver,SolvType);
Prob            = iniSolve(Prob,SolvType,1,1);
n               = Prob.N;

PriLev          = Prob.PriLev;             % Print level
fGoal           = Prob.optParam.fGoal;     % Goal for f(x).
fTol            = Prob.optParam.eps_f;     % Relative tolerance for fGoal

if isempty(fGoal), fGoal = -inf; end

x_L             = Prob.x_L(:);             % Lower bounds on x
x_U             = Prob.x_U(:);             % Upper bounds on x
bTol            = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol            = Prob.optParam.cTol;      % Nonlinear constraint feasibility tolerance
MaxFunc         = Prob.optParam.MaxFunc;   % Max # of function evaluations in each glcDirect subproblem
MaxIter         = Prob.optParam.MaxIter;   % Max # of iterations in each glcDirect subproblem
IntVars         = DefPar(Prob.MIP,'IntVars',[]);  % Integer variables
IV              = false(n,1);              % Logical vector for integers
dLin            = Prob.mLin;
dNoL            = Prob.mNonLin;
K               = dLin+dNoL;

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('multiMINLP: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars         = find(IV);
RealVars        = find(~IV);
eps_I           = 1E-8;
nI              = length(IntVars);
nR              = length(RealVars);
xD              = x_U-x_L;
if nI > 0 
   nMax  = DefPar(Prob.MIP,'nMax',[]);
   if isempty(nMax) | ~isfinite(nMax)
      nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
   end
else
   nMax  = 0;
end
MaxFunc               = max(MaxFunc,MaxIter);
MaxIter               = MaxFunc;
Prob.optParam.MaxFunc = MaxFunc;
Prob.optParam.MaxIter = MaxIter;

xInit = DefPar(Prob,'xInit',[]);

if nI > 0
   if isempty(xInit), xInit = min(3000,30*n); end
   if nMax < 300 & xInit < nMax
      xInit = nMax*2;
   end
end

if nI == 0 | (nI > 0 & nI < n & nMax <= xInit)
   % fprintf('multiMINLP: No integer variables. Call multiMin instead\n');
   Prob.xInit             = xInit;
   Result                 = tomRun('multiMin',Prob, PriLev-4);
   Result.Solver = ['multiMINLP calling ',Result.Solver];
   Result.multiMINLP.nMax = 0;
   Result.multiMINLP.Type = 4;
   return
end

Rfac  = DefPar(Prob.MIP,'Rfac',0.25);

if Prob.WarmStart == 1 & isfield(Prob,'multiMINLP')
   MINLP  = Prob.multiMINLP;
   Type   = DefPar(MINLP,'Type',0);
   if Type == 1
      WarmStart   = 1;
      MM          = Prob.multiMin;
      nMax        = MINLP.nMax;
      nComp       = MINLP.nComp;
      nTry        = MINLP.nTry;
      Tmax        = MINLP.Tmax;
      T           = MINLP.T;
      % Should classify if the T:th region fully searched or not
      % Possibly shrink the region to the one actually searched
      XL          = MINLP.XL(:,1:T);
      XU          = MINLP.XU(:,1:T);
      Runs        = MINLP.Runs;
      Feasible    = MINLP.Feasible;
      fIP         = MINLP.fIP;
      fPen        = MINLP.fPen;
      xIP         = MINLP.xIP;
      % HKH - for now set WarmStart == 0 
      WarmStart   = 0;
   else
      WarmStart   = 0;
   end
   Prob.WarmStart = 0;  % Must set 0 to avoid warm start in NLP solver
else
   WarmStart      = 0;
   Prob.WarmStart = 0;
end


if (nMax <= MaxFunc & nI == n ) % | ( nMax <= min(MaxFunc,5000) & nI <= n-2 )
   % Pure IP with few (<= MaxFunc) IP combinations
   Result   = tomRun('glcDirect',Prob,PriLev-4);
   Result.Solver           = ['multiMINLP calling glcDirect'];
   Result.SolverAlgorithm  = 'Find IP optimum with multiMINLP algoritm';
   Result.multiMINLP.Type  = 3;
   Result.multiMINLP.nMax  = nMax;
   Result.multiMINLP.nComp = Result.glcDirect.WarmStartInfo.nFuncTot;
   Result.multiMINLP.nTry  = Result.glcDirect.WarmStartInfo.nRectTot;
   Result.multiMin         = [];
   Result                  = endSolve(Prob,Result);
   Result.f_k              = Result.f_k-Prob.fConstant;
   return
else
   ProbMM                  = Prob;
   ProbMM.MIP.IntVars      = [];
   rMM                     = tomRun('multiMin',ProbMM,PriLev-4);
   if rMM.ExitFlag == 4
    Result                 = rMM;
    Result.Solver          = ['multiMINLP - call multiMin, Local solver ' localSolver];
    Result.multiMINLP.nMax = nMax;
    Result.multiMINLP.Type = 2;
   else
    Info     = rMM.multiMin.Info;
    FuncEv   = Info.FuncEv;
    GradEv   = Info.GradEv;
    HessEv   = Info.HessEv;
    ConstrEv = Info.ConstrEv;
    nGlobal  = Info.nGlobal;
    nLocal   = Info.nLocal;
    fOpt     = rMM.multiMin.FX(1,1:nLocal);
    xOpt     = rMM.multiMin.XX(:,1:nLocal);
    ixOpt    = 1:nLocal;
    Icomp    = n*ones(1,nLocal);
    X        = rMM.x_k;
    % Step 1: Pick best minimum and solve around this minimum
    [IC,Iidx]       = IntComp(xOpt(:,1:nGlobal),IntVars,eps_I);
    if length(IC) > 1
       % Select the x with max integer components
       [icMax, idx] = max(IC);
       x_T          = xOpt(:,idx);
       Ridx         = IntVars(~Iidx(:,idx));
       Icomp(idx)   = icMax;
    else
       idx          = 1;
       Icomp(idx)   = IC;
       Ridx         = IntVars(~Iidx);
    end
    ProbIP          = Prob;
    x_T             = xOpt(:,idx);

    % Restrict real variables
    if nR > 0
       ProbIP.x_L(RealVars) = max(x_T(RealVars)-Rfac*abs(x_T(RealVars)),x_L(RealVars));
       ProbIP.x_U(RealVars) = min(x_T(RealVars)+Rfac*abs(x_T(RealVars)),x_U(RealVars));
    end
    % Fix integer valued integer variables
    Ii              = IntVars(Iidx(:,idx));
    x_T(Ii)         = round(x_T(Ii));
    ProbIP.x_L(Ii)  = x_T(Ii);
    ProbIP.x_U(Ii)  = x_T(Ii);
    % How many integer combinations left among non-integer valued integer variables
    nSMax           = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
    nSMaxM          = nSMax;
    if nSMax <= MaxIter
       % Keep the bounds on the real valued integer variables
       s               = inf;
    else
       % Reduce the number of combinations
       s                = max(1,floor((exp(log(MaxIter)/length(Ridx))-2)/2));
       % Restrict integer variables
       ProbIP.x_L(Ridx) = max(floor(x_T(Ridx))-s,x_L(Ridx));
       ProbIP.x_U(Ridx) = min(ceil(x_T(Ridx))+s,x_U(Ridx));
       nSMax           = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
       while s > 0 & nSMax > MaxIter
           s = s-1;
           ProbIP.x_L(Ridx) = max(floor(x_T(Ridx))-s,x_L(Ridx));
           ProbIP.x_U(Ridx) = min(ceil(x_T(Ridx))+s,x_U(Ridx));
           nSMax           = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
       end
       if nSMax > MaxIter
          x_I              = floor(x_T(Ridx)+eps_I*max(1,abs(x_T(Ridx))));
          x_r              = max(0,x_T(Ridx)-x_I);
          r                = abs(x_r-0.5);
          [rSmall iSmall]  = sort(r,'ascend');
          j                = 0;
          while nSMax > MaxIter
             j                = j+1;
             idx              = Ridx(iSmall(j));
             iVal             = ProbIP.x_U(idx)-ProbIP.x_L(idx)+1;
             ProbIP.x_L(idx) = round(x_T(idx));
             ProbIP.x_U(idx) = round(x_T(idx));
             nSMax            = nSMax/iVal;
          end
       end
       nSMax = prod(1+(ProbIP.x_U(IntVars)-ProbIP.x_L(IntVars)));
    end

    ProbIP.x_0 = x_T;
    if nSMax <= 5*MaxIter
       rIP                       = tomRun('glcFast',ProbIP,PriLev-4);
       % rIP                     = tomRun('glcSolve',ProbIP,PriLev-4);
       % rIP                     = tomRun('glcDirect',ProbIP,PriLev-4);
    else 
       ProbIP.optParam.IterPrint = 1;
       rIP                       = tomRun('minlpSolve',ProbIP,PriLev-4);
    end
    [hL1 h]  = consviolation(rIP, PriLev-5, 0);
    Inform   = rIP.Inform;
    ExitFlag = rIP.ExitFlag;
%   if (hL1 >= K*cTol & Inform == 4) | (Inform ~= 92 & nR == 0  & nTry < nSMax)
%      ProbIP.WarmStart = 1;
%      rIP              = tomRun('glcFast',ProbIP,PriLev-4);
%      Inform           = rIP.Inform;
%      ExitFlag         = rIP.ExitFlag;
%      %ProbIP          = WarmDefGLOBAL('glcDirect',ProbIP,rIP);
%      %rIP             = tomRun('glcDirect',ProbIP,PriLev-4);
%      [hL1 h]          = consviolation(rIP, PriLev-5, 0);
%   end
%   ProbIP.WarmStart = 0;

    FuncEv    = FuncEv   + rIP.FuncEv;
    GradEv    = GradEv   + rIP.GradEv;
    HessEv    = HessEv   + rIP.HessEv;
    ConstrEv  = ConstrEv + rIP.ConstrEv;
    nComp     = rIP.FuncEv;
    nTry      = rIP.Iter;
    Feasible  = hL1 <= 100*K*cTol;
    fPen      = rIP.f_k+hL1;
    if Feasible
       fIP    = fPen;
       xIP    = rIP.x_k(:,1);
    else
       fIP    = inf;
       xIP    = [];
    end
    Runs      = idx;
    ixOpt     = setdiff(ixOpt,idx);
    % Next Steps
    if ~isinf(fIP)
       ixOpt  = setdiff(ixOpt,ixOpt( fIP*(1+1E-12) < fOpt(ixOpt) ));
    end
    Tmax      = 1+length(ixOpt);
    if PriLev > 0
       fprintf('   TotIP %d. IP comb %d. Reduced to %d. s=%d.\n',nMax,nSMaxM,nSMax,s);
       fprintf('#  1-');
       fprintf('%3d ',idx);
       fprintf('fIP %12.7f. ',fIP);
       fprintf('fPen %12.7f. ',fPen);
       fprintf('fOpt %12.7f. ',fOpt(idx));
       fprintf('f_k %12.7f. ',rIP.f_k);
       fprintf('hL1%13.7f. ',hL1);
       fprintf('I%3d. ',Icomp(idx));
       fprintf('R%3d. ',length(Ridx));
       fprintf('Exit %2d. ',ExitFlag);
       fprintf('Inform %2d. ',Inform);
       fprintf('It%5d. ',rIP.Iter);
       fprintf('FuEv%5d. ',rIP.FuncEv);
       fprintf('Comp f(x)%6d. ',nComp);
       fprintf('Check %d other min ',Tmax-1);
       fprintf('\n');
    end
    XL        = zeros(n,Tmax);
    XU        = zeros(n,Tmax);
    T         = 1;
    XL(:,1)   = ProbIP.x_L;
    XU(:,1)   = ProbIP.x_U;
    if Tmax > 1
       % Search around other local minima that has fOpt < fIP
       while ~isempty(ixOpt)
          k                   = ixOpt(1);
          inside              = 0;
          x_T                 = xOpt(:,k);
          [IC,Iidx]           = IntComp(x_T,IntVars,eps_I);
          Icomp(k)            = IC;
          Ridx                = IntVars(~Iidx);
          % Round integer valued integer variables
          Ii                  = IntVars(Iidx);
          x_T(Ii)             = round(x_T(Ii));
          for j=1:T                  % Check if new point in area already searched
              if all(XL(:,j) <= x_T) & all(XU(:,j) >= x_T)
                 inside       = 1;
                 ixOpt        = setdiff(ixOpt,k);
                 break;
              end
          end
          if ~inside
             ProbIP.x_L       = x_L;
             ProbIP.x_U       = x_U;
             % Restrict real variables
             if nR > 0
                ProbIP.x_L(RealVars) = max(x_T(RealVars)-Rfac*xD(RealVars),x_L(RealVars));
                ProbIP.x_U(RealVars) = min(x_T(RealVars)+Rfac*xD(RealVars),x_U(RealVars));
             end
             ProbIP.x_L(Ii)   = x_T(Ii);
             ProbIP.x_U(Ii)   = x_T(Ii);
             % How many integer combinations left
             nSMax            = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
             nSMaxM           = nSMax;
             T                = T+1;
             if nSMax <= MaxIter
                % Keep the bounds on the real valued variables
                s                = inf;
             else
                % Reduce the number of combinations
                s                = max(1,floor((exp(log(MaxIter)/length(Ridx))-2)/2));
                % Restrict integer variables
                ProbIP.x_L(Ridx) = max(floor(x_T(Ridx))-s,x_L(Ridx));
                ProbIP.x_U(Ridx) = min(ceil(x_T(Ridx))+s,x_U(Ridx));
                nSMax           = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
                while s > 0 & nSMax > MaxIter
                   s = s-1;
                   ProbIP.x_L(Ridx) = max(floor(x_T(Ridx))-s,x_L(Ridx));
                   ProbIP.x_U(Ridx) = min(ceil(x_T(Ridx))+s,x_U(Ridx));
                   nSMax            = prod(1+(ProbIP.x_U(Ridx)-ProbIP.x_L(Ridx)));
                end
                if nSMax > MaxIter
                   x_I              = floor(x_T(Ridx)+eps_I*max(1,abs(x_T(Ridx))));
                   x_r              = max(0,x_T(Ridx)-x_I);
                   r                = abs(x_r-0.5);
                   [rSmall iSmall]  = sort(r,'ascend');
                   j                = 0;
                   while nSMax > MaxIter
                      j                = j+1;
                      idx              = Ridx(iSmall(j));
                      iVal             = ProbIP.x_U(idx)-ProbIP.x_L(idx)+1;
                      ProbIP.x_L(idx)  = round(x_T(idx));
                      ProbIP.x_U(idx)  = round(x_T(idx));
                      nSMax            = nSMax/iVal;
                   end
                end
                nSMax               = prod(1+(ProbIP.x_U(IntVars)-ProbIP.x_L(IntVars)));
             end
             if ~isinf(fIP)
                ProbIP.MIP.fIP      = fIP;
                ProbIP.MIP.xIP      = xIP;
             end
             if nSMax <= 5*MaxIter
                % Find best integer point close to minimum
                rIPi                      = tomRun('glcFast',ProbIP,PriLev-4);
                %rIPi                     = tomRun('glcDirect',ProbIP,PriLev-4);
             else 
                ProbIP.optParam.IterPrint = 1;
                rIPi                      = tomRun('minlpSolve',ProbIP,PriLev-4);
             end

             FuncEv              = FuncEv   + rIPi.FuncEv;
             GradEv              = GradEv   + rIPi.GradEv;
             HessEv              = HessEv   + rIPi.HessEv;
             ConstrEv            = ConstrEv + rIPi.ConstrEv;
             nComp               = nComp    + rIPi.FuncEv;
             nTry                = nTry     + rIPi.Iter;
             [hL1 h]             = consviolation(rIPi, PriLev-5, 0);
             fPeni               = rIPi.f_k + hL1;
             if hL1 <= 100*K*cTol
                fIPi             = fPeni;
             else
                fIPi             = inf;
             end
             Runs                = [Runs,k];
             ixOpt               = setdiff(ixOpt,k);
             if fIPi < fIP
                % Update of best IP point found
                if PriLev > 1
                   fprintf('      ');
                   fprintf('Local IP min around local k = %d better \n',k);
                end
                Feasible = 1;
                if PriLev > 1
                   fprintf('      ');
                   fprintf('Old fIP %12.7f. New fIP %12.7f\n',fIP,fIPi);
                end
                fIP      = fIPi;
                xIP      = rIPi.x_k(:,1);
                ixOpt    = setdiff(ixOpt,ixOpt( fIP*(1+1E-12) < fOpt(ixOpt) ));
                rIP      =    rIPi;
             elseif ~Feasible & fPeni < fPen
                if PriLev > 1
                   fprintf('      ');
                   fprintf('Old fPen %12.7f. New fPen %12.7f\n',fPen,fPeni);
                end
                fPen     = fPeni;
                xIP      = rIPi.x_k(:,1);
                rIP      = rIPi;
                if ~isinf(fIP), keyboard; end
             end
             % Update with new region searched through
             XL(:,T) = ProbIP.x_L;
             XU(:,T) = ProbIP.x_U;
             if PriLev > 0
                fprintf('   TotIP %d. IP comb %d. Reduced to %d. s=%d.\n',nMax,nSMaxM,nSMax,s);
                fprintf('#%3d-',T);
                fprintf('%3d ',k);
                fprintf('fIP %12.7f. ',fIP);
                fprintf('fPen %12.7f. ',fPen);
                fprintf('fOpt %12.7f. ',fOpt(k));
                fprintf('f_k %12.7f. ',rIPi.f_k);
                fprintf('hL1%13.7f. ',hL1);
                fprintf('I%3d. ',Icomp(k));
                fprintf('R%3d. ',length(Ridx));
                fprintf('Exit %2d. ',ExitFlag);
                fprintf('Inform %2d. ',Inform);
                fprintf('It%5d. ',rIPi.Iter);
                fprintf('FuEv%5d. ',rIPi.FuncEv);
                fprintf('Comp f(x)%6d. ',nComp);
                fprintf('\n');
             end
          end
       end
       if PriLev > 0
          fprintf('    ');
          fprintf('fIP %12.7f. ',fIP);
          fprintf('Solved %d IP problems out of %d. ',T,length(fOpt));
          fprintf('FuncEv %d GradEv %d ConstrEv %d',FuncEv, GradEv,ConstrEv);
          fprintf('\n');
          fprintf('    ');
          fprintf('Total IP-comb: %20.0f. ',nMax);
          fprintf('Comp f(x) %d. ',nComp);
          fprintf('\n');
       end
    end
    Result                     = rIP;
    Result.multiMINLP.Type     = 1;
    Result.multiMINLP.nMax     = nMax;
    Result.multiMINLP.nComp    = nComp;
    Result.multiMINLP.nTry     = nTry;
    Result.multiMINLP.Tmax     = Tmax;
    Result.multiMINLP.T        = T;
    Result.multiMINLP.XL       = XL(:,1:T);
    Result.multiMINLP.XU       = XU(:,1:T);
    Result.multiMINLP.Runs     = Runs;
    Result.multiMINLP.Feasible = Feasible;
    Result.multiMINLP.fIP      = fIP;
    Result.multiMINLP.fPen     = fPen;
    Result.multiMINLP.xIP      = xIP;
    Result.multiMin            = rMM.multiMin;
    Result.FuncEv              = FuncEv;
    Result.GradEv              = GradEv;
    Result.HessEv              = HessEv;
    Result.ConstrEv            = ConstrEv;
    Result.Solver              = ['multiMINLP - glcDirect and multiMin, Local solver ' localSolver];
   end
   Result.SolverAlgorithm      = 'multiMINLP Hybrid Algoritm';
   Result                      = endSolve(Prob,Result);
   Result.f_k                  = Result.f_k-Prob.fConstant;
   return
end

% NOT USED NOW
function ix = checkEqual(x_k,xOpt,xEqTol)
ix = 0;
for i=1:size(xOpt,2)
    if norm(x_k-xOpt(:,i)) <= xEqTol*max(1,norm(x_k))
        ix = i;
        break;
    end
end

% WRONG?
function h_k = consViol(Prob,Result,bTol,cTol)
x_k = Result.x_k; 
h_k = 0;
if Prob.mLin > 0
   L = Prob.A*x_k;
   h_k = h_k+sum(max(0,max(Prob.b_L-bTol*max(1,abs(Prob.b_L))-L,...
                      L-bTol*max(1,abs(Prob.b_U))-Prob.b_U)));
end
if Prob.mNonLin > 0
   C = Result.c_k;
   h_k = h_k+sum(max(0,max(Prob.c_L-cTol*max(1,abs(Prob.c_L))-C,...
                      C-cTol*max(1,abs(Prob.c_U))-Prob.c_U)));
end

% MODIFICATION LOG:
%
% 090911  hkh  Written
% 091002  hkh  Use WarmDefGLOBAL to warm start glcDirect
% 091017  hkh  Revision due to revised multiMin
% 091021  hkh  Prepared code for warm start, but no algorithm yet for this
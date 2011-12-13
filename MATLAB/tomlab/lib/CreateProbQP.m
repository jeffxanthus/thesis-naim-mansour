% CreateProbQP creates a basic Prob structure for a Quadratic Programming (QP)
% subproblem, as well as LP, Feasible Point (FP), and Dual LP (DLP).
%
% -----------------------------------------------------
%
% QP minimization problem:
%
%
%        min   0.5 * x' * F * x + c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of CreateProbQP:
%
% function ProbQP = CreateProbQP(Prob, Type, MaxIter, PriLev, wait)
%
% INPUT (One parameter F must always be given. Empty gives default)
%
% Prob       Input problem structure
%            If Prob.QP.SOL or Prob.QP.optParam are defined, then
%            these structures are moved to ProbQP.SOL and ProbQP.optParam.
%
% Type       Any of QP, LP, FP, DLP
% MaxIter    Maximal number of iterations
% PriLev     Print level
% wait       Flag if to pause each iteration (only in m-file solvers)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Nov 6, 2000. Last modified Jul 17, 2009.

function ProbQP = CreateProbQP(Prob, Type, MaxIter, PriLev, wait)

if nargin < 5
   wait=[];
   if nargin < 4
      PriLev=[];
      if nargin < 3
         MaxIter=[];
         if nargin < 1
            error('CreateProbQP requires two parameters Prob and Type');
         end
      end
   end
end

ProbQP = ProbDef(0);
Type   = lower(Type);

switch Type(1:1)
   case 'q'
      probType        = checkType('qp');
      Solver          = Prob.SolverQP;
      ProbQP.Name     = [Prob.Name,' QP subproblem'];
   case 'l'
      probType        = checkType('lp');
      Solver          = Prob.SolverLP;
      ProbQP.Name     = [Prob.Name,' LP subproblem'];
   case 'f'
      probType        = checkType('lp');
      Solver          = Prob.SolverFP;
      ProbQP.Name     = [Prob.Name,' FP subproblem'];
   case 'd'
      probType        = checkType('lp');
      Solver          = Prob.SolverDLP;
      ProbQP.Name     = [Prob.Name,' DLP subproblem'];
end

n           = Prob.N;
ProbQP.N    = n;
ProbQP.P    = Prob.P;
ProbQP.x_0  = zeros(n,1);
ProbQP.x_L  = -Inf*ones(n,1);
ProbQP.x_U  = Inf*ones(n,1);
ProbQP.A    = Prob.A;
ProbQP.b_L  = Prob.b_L;
ProbQP.b_U  = Prob.b_U;
ProbQP.mLin = size(Prob.A,1);
ProbQP.mNonLin = 0;

ProbQP.PriLevOpt        = PriLev;
ProbQP.optParam.MaxIter = MaxIter;
ProbQP.optParam.wait    = wait;
ProbQP.probFile         = 0;
ProbQP.probType         = probType;
ProbQP.SolverQP         = Prob.SolverQP;
ProbQP.SolverLP         = Prob.SolverLP;
ProbQP.SolverDLP        = Prob.SolverDLP;
ProbQP.SolverFP         = Prob.SolverFP;
ProbQP.f_Low            = Prob.f_Low;
ProbQP.LargeScale       = Prob.LargeScale;
ProbQP.solvType         = checkType('qp');
ProbQP.CHECK            = 1;

if isfield(Prob.QP,'Solver')
   ProbQP.Solver=Prob.QP.Solver;
else
   ProbQP.Solver.Alg=0;
   ProbQP.Solver.Method=0;
   ProbQP.Solver.Name=Solver;
end

if ~isfield(Prob.QP,'SOL')
    ProbQP.SOL = struct( 'SpecsFile', [], 'PrintFile', [], 'SummFile', [], ...
        'xs', [], 'hs', [], 'nS', 0, 'hElastic', [], ...
        'iState', [], 'cLamda', [], 'R', [], 'optPar', -999*ones(1,71), ...
        'optParN', 71);
else
    ProbQP.SOL = Prob.QP.SOL;
end

if ~isfield(Prob.QP,'optParam')
   ProbQP.optParam = optParamDef(Solver,probType, n, 0, size(Prob.A,1));
else
   % Check the optParam structure
   ProbQP.optParam = optParamSet(Prob.QP.optParam,Solver,probType);
end
% Increase number of iterations, sometimes the default is too low
ProbQP.optParam.MaxIter = max(ProbQP.optParam.MaxIter, MaxIter);

ProbQP.QP = struct( 'F', [], 'c', [], 'B', [], 'y', [], ...
                    'Q', [], 'R', [], 'E', [], ...
                    'Ascale', [], 'DualLimit', [], ...
                    'UseHot', [], 'HotFile', [], 'HotFreq', [], 'HotN', []);

if Type(1:1) == 'q'
    ProbQP = tomFiles(ProbQP, 'qp_f', 'qp_g', 'qp_H');
else
    ProbQP = tomFiles(ProbQP, 'lp_f', 'lp_g', 'lp_H');
end

% MODIFICATION LOG
%
% 001106  hkh  Written
% 010903  hkh  Change to 63 optPar parameters
% 020304  hkh  Increase maximal iterations to max(MaxIter,default).
% 020524  hkh  Call mFiles, not set FUNCS structure explicit
% 040111  hkh  Define fields mLin and mNonLin
% 040728  med  tomFiles used instead
% 050605  hkh  Change to optParN = 65, for SNOPT 6
% 050606  hkh  Change to optParN = 71, for SNOPT 7
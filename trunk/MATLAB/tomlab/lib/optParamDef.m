% Initialization of structure optParam
%
% The fields in optParam store different optimization parameters
% that control the operation of the solvers, e.g. convergence tests
%
%
% function [optParam, SOLPar] = optParamDef(Solver,probType,nObj,nJac,m);
%
% INPUT:
%   Solver   Name of solver
%   probType Problem type
%   nObj     Number of nonlinear objective variables
%   nJac     Number of nonlinear constraint variables
%   m        Number of rows in the constraint matrix (linear & nonlinear)
%
% OUTPUT:
%   optParam Structure
%   SOLPar   The SOL parameter vector if a SOL server, otherwise empty

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written Sept 11, 1998.  Last modified March 10, 2008.

function [optParam, SOLPar] = optParamDef(Solver,probType,nObj,nJac,m)

% Base some of the defaults on the defaults in SNOPT (similar to MINOS)
if nargin < 5
   m = [];
   if nargin < 4
      nJac = [];
      if nargin < 3
         nObj = [];
         if nargin < 2
            probType = 3;
            if nargin < 1
               Solver = 'snopt';
            end
         end
      end
   end
end

MaxFunc   = 10000;
IterPrint = 0;
EpsF      = 1E-8;
SOLSolver = 0;
WAIT      = 0;
PRESOLVE  = 0;
xTol      = 1000*eps;
if findstr(Solver,'Solve')
   switch Solver(1:3)
     case 'ucS'
       SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
       SOLPar(30)=500;    % MaxIter
     case 'qpS'
       SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
       SOLPar(30)=2000;   % MaxIter
       SOLPar(11)=1E-8;   % bTol
     case 'con'
       SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
       SOLPar(30)=500;    % MaxIter
       SOLPar(11)=1E-8;   % bTol
     case 'nlp'
       SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
       SOLPar(30)=500;    % MaxIter
       SOLPar(11)=1E-8;   % bTol
     case 'mip'
       SOLPar = SOLGet('lpopt',probType,nObj,nJac,m);
       SOLPar(30)=5000;   % MaxIter
       SOLPar(11)=1E-8;   % bTol
       IterPrint = 1;
       SOLPar(1) = 0;
     case 'minlp'
       % New solver minlpSolve
       SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
       SOLPar(30)=2000;    % MaxIter
       SOLPar(11)=1E-8;    % bTol
     case 'cls'
       SOLPar = SOLGet('nlssol',probType,nObj,nJac,m);
       SOLPar(30)=500;          % MaxIter
       SOLPar(11)=1E-8;         % bTol
       SOLPar(27)=eps^(0.67);   % eps_Rank
       SOLPar(42)=3.0E-13^(0.5);% snopt and minos DiffInt 
       % Central differences, SOLPar(43) is also not set, but not used 
     case 'lpS'
       SOLPar = SOLGet('lpopt',probType,nObj,nJac,m);
       SOLPar(30)=2000;   % MaxIter
       SOLPar(11)=1E-8;   % bTol
     case 'glb'
       SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
       % MaxFunc    = 1000;
       MaxFunc    = max(10000,nObj*2000);
       SOLPar(30) = MaxFunc;
       EpsF       = 1E-2;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
     case 'glc'
       SOLPar = SOLGet('npsol',probType,nObj,nJac,m);
       SOLPar(30) = 1000;   % MaxIter
       % MaxFunc    = 1000;
       MaxFunc    = max(10000,nObj*2000);
       SOLPar(9)  = 1E-5;   % cTol
       SOLPar(11) = 1E-7;   % bTol
       SOLPar(10) = 1E-11;  % eps_x
       EpsF       = 1E-2;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
     case 'rbf'
       SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
       SOLPar(30) = 10000;
       EpsF       = 1E-3;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
       SOLPar(9)  = 1E-5;   % cTol
       SOLPar(11) = 1E-7;   % bTol
       MaxFunc    = 200;
       IterPrint  = 1;
     otherwise
       SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   end
elseif strcmpi(Solver,'arbfmip') | strcmpi(Solver,'arbf')
   SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
   SOLPar(27) = eps^(0.67);   % eps_Rank
   SOLPar(30) = 10000;    % Max iterations
   EpsF       = 1E-3;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
   SOLPar(9)  = 1E-5;   % cTol
   SOLPar(11) = 1E-7;   % bTol
   MaxFunc    = 200;
   IterPrint  = 1;
elseif strcmpi(Solver,'ego')
   %SOLPar(1:63) = -999;
   SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
   SOLPar(27) = eps^(0.67);   % eps_Rank
   SOLPar(30) = 10000;    % Max iterations
   SOLPar(9)  = 1E-5;     % cTol
   SOLPar(11) = 1E-7;     % bTol
   MaxFunc    = 200;      % Restrict only on function values, not iterations
   EpsF       = 1E-2;     % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
   IterPrint  = 1;
elseif strcmpi(Solver,'cutplane')
   SOLPar = SOLGet('lpopt',probType,nObj,nJac,m);
   SOLPar(1) = 0;
   SOLPar(30)=200;
   IterPrint = 1;
elseif strcmpi(Solver,'sTrustr')
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=500;
elseif strcmpi(Solver,'qld')
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
elseif strcmpi(Solver,'bqpd')
   % Dundee solver 
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'Xpress-MP')
   % Xpress-MP from Dash Optimization
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'CPLEX')
   % CPLEX from ILOG
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'XA')
   % XA 
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'MIQPbb')
   % Dundee solver 
   SOLPar = SOLGet('qpopt',probType,nObj,nJac,m);
   SOLPar(30)= 20000;
   SOLPar(10)= 1E-10;
   EpsF      = 1E-5;
elseif strcmpi(Solver,'filterSQP')
   % Dundee solver 
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'nlpql')
   % KS solver 
   SOLPar = SOLGet('npsol',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'PDCO')
   % SOL convex NLP solver 
   options = pdcoSet;
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=options.MaxIter;       % optParam.Maxiter
   SOLPar(36)=options.LSQRMaxIter;   % optParam.MinorIter
   SOLPar(10)=options.OptTol;        % optParam.eps_x
   SOLPar(11)=options.FeaTol;        % optParam.bTol
   WAIT      =options.wait;          % optParam.wait
elseif strcmpi(Solver,'PDSCO')
   % SOL convex NLP solver 
   options = pdcoSet;
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=options.MaxIter;       % optParam.Maxiter
   SOLPar(36)=options.LSQRMaxIter;   % optParam.MinorIter
   SOLPar(10)=options.OptTol;        % optParam.eps_x
   SOLPar(11)=options.FeaTol;        % optParam.bTol
   WAIT      =options.wait;          % optParam.wait
elseif strcmpi(Solver,'MINLPbb')
   % Dundee solver 
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=5000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'PENSDP')
   % PENSDP
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'PENBMI')
   % PENBMI
   SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
   SOLPar(30)=2000;
   SOLPar(10)=1E-10;
elseif strcmpi(Solver,'glbFast') | strcmpi(Solver, 'glbDirect')
   SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
   % SOLPar(30) = 200;
   % MaxFunc    = 10000;
   MaxFunc    = max(10000,nObj*2000);
   SOLPar(30) = MaxFunc;
   EpsF       = 1E-2;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
elseif strcmpi(Solver,'glcFast') | strcmpi(Solver, 'glcDirect')
   SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
   % MaxFunc    = 10000;
   MaxFunc    = max(10000,nObj*2000);
   SOLPar(30) = MaxFunc;
   % SOLPar(30) = 10000;
   EpsF       = 1E-2;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
   SOLPar(9)  = 1E-5;   % cTol
   SOLPar(11) = 1E-7;   % bTol
   SOLPar(10) = 1E-11;  % eps_x
elseif strcmpi(Solver,'glcCluster')
   SOLPar     = SOLGet('npsol',probType,nObj,nJac,m);
   SOLPar(30) = 200;
   EpsF       = 1E-2;   % Relative tolerance (f(x) - fOpt) <= EpsF * fOpt
   SOLPar(11) = 1E-7;   % bTol
   SOLPar(9)  = 1E-5;   % cTol
   SOLPar(10) = 1E-11;  % eps_x
   %MaxFunc    = 200;
   MaxFunc    = 10000;
elseif strcmpi(Solver,'lsei')

   SOLPar = SOLGet('lssol',probType,nObj,nJac,m);
   SOLPar(30)=500;
   SOLPar(1) = 0;
elseif strcmpi(Solver,'constr')

   SOLPar = SOLGet('npsol',probType,nObj,nJac,m);
   SOLPar(30)=1000;
   SOLPar(1) = 0;
elseif strcmpi(Solver,'Tlsqr')
   SOLPar = SOLGet('lssol',probType,nObj,nJac,m);
   SOLPar(30)=1000;
   SOLPar(1) = 0;
else
   % A SOL solver
   SOLPar = SOLGet(Solver,probType,nObj,nJac,m);
   if isempty(SOLPar)
      % An opt tbx solver
      SOLPar = SOLGet('snopt',probType,nObj,nJac,m);
      SOLPar(30)=10000;
   else
      SOLSolver = 1;
   end
   if strcmpi(Solver,'minos')
      % Will influence printing in PrintResult
      xTol = 1E-10;
      %alfa=SOLPar(22)
   end
end

optParam = struct( ...
  'PriLev',SOLPar(1), 'PriFreq',SOLPar(5), 'SummFreq',SOLPar(6), ...
  'MinorPriLev',SOLPar(2), ...
  'IterPrint',IterPrint, 'wait',WAIT, 'MaxFunc',MaxFunc, ... 
  'MaxIter',SOLPar(30), 'MajorIter',SOLPar(35), 'MinorIter',SOLPar(36), ...
  'eps_f',EpsF, 'eps_absf',1000*eps, 'eps_x',SOLPar(10), 'eps_dirg',1E-8, ...
  'eps_g',1E-7, 'eps_Rank',SOLPar(27), 'EpsGlob', 1E-4,'fGoal',-Inf,...
  'fTol',SOLPar(41), 'xTol',xTol, 'bTol',SOLPar(11), 'cTol',SOLPar(9),...
  'MinorTolX',SOLPar(12), ...
  'size_x',1, 'size_f',1, 'size_c',1, 'PreSolve', PRESOLVE, ...
  'DerLevel',SOLPar(39), 'GradCheck',  SOLPar(13), ...
  'DiffInt', 1E-6, 'CentralDiff',1E-4, ...
  'QN_InitMatrix',[], ...
  'splineSmooth', -1, 'splineTol', 1E-3, ...
  'BigStep',SOLPar(45), 'BigObj',SOLPar(46), 'CHECK',0);

if SOLPar(42) > 0
   optParam.DiffInt = SOLPar(42);
end
if SOLPar(43) > 0
   optParam.CentralDiff = SOLPar(43);
end

if nargout > 1
   if ~SOLSolver
      SOLPar = -900*ones(1,72);
   end
end

return

% Some comments about the optimization parameters

% (x) refers to the OPT TBX 1.x options vector

% (2) eps_x: Convergence tolerance in optimal solution x
% Relative distance between successive x. //x_k+1 - x_k//. Def OPTIM 1E-4.

% (3) eps_dirg: Convergence tolerance on F. Def OPTIM 1E-4. Here 1E-8.
% Directed.derivative:  g_k^T * p_k <= eps_dirg 

% (4) cTol (eps_c): Constraint violation convergence tolerance. Def OPTIM 1E-6.

% OPTIONS(8)  - Function value. (Lambda in goal attainment. )
% Returned as Result.f_k

% (9) GracCheck: Check user-supplied gradients, -1 no check, 
%     0-3 different VERIFY LEVELS

% OPTIONS(10) - Number of Function and Constraint Evaluations.
% OPTIONS(11) - Number of Function Gradient Evaluations.
% OPTIONS(12) - Number of Constraint Evaluations
% OPTIONS(13) - Number of equality constraints. 
% This number is implicit, i.e.
% me = OPTIONS(13)= sum(Prob.b_L==Prob.b_U) +  sum(Prob.c_L==Prob.c_U)

% (14) MaxIter: Maximal number of iterations

% OPTIONS(15) - Used in goal attainment for special objectives. 

% (16) Minimum change in variables for finite difference gradients.
% Def OPTIM 1E-8. Here 1E-8.

% (17) Maximum change in variables for finite difference gradients.
% Def OPTIM 0.1. Here 0.1.


% (20) Gradient (or reduced gradient) convergence tolerance
%  max_i | g_k(i)*max(x(i),size_x) | < eps_g * max( |f(x)|, size_f)
% Squared before calling SOL, as they take square root first.


% (22) Lower bound on function value. Used in line search by Fletcher.
% Set at top level of struct Prob: Prob.f_Low=-1E300;

% (23) eps_Rank: Rank test tolerance. PIVOT TOLERANCE in SOL.

% (24) wait: Flag if to use pause statements after output

% (25) Convergence tolerance on Absolute function value     |f_k| < eps_absf

%     OPTIONS(26) - iter k. Number of main (major) iterations.
%     OPTIONS(27) - Number of minor iterations.
%     OPTIONS(28) - EXIT flag, convergence to local min = 0. Otherwise errors.
%     OPTIONS(29) - INFORM, information parameter, type of convergence.

% New parameters in optParam

% PreSolve: Flag if presolve analysis is to be applied on linear constraints

% Initial matrix for Quasi-Newton, may be set by the user.
% When QN_InitMatrix is empty, the identity matrix is used.

% xTol: If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bound
% bTol: Linear feasibility tolerance
% optParam.bTol=1E-8;    
% cTol: Nonlinear feasibility tolerance
% optParam.cTol=1E-6;    
% fTol: Accuracy in the computation of the function value. 
% Assume eps^0.8 (Same as SOL routines - FUNCTION PRECISION)

% size_x, size_f, size_c:
% Size at optimum for variables x, function f and constraints c
% Used to get good relative convergence tests. Default 1.

% eps_f: If (f_km1-f_k) < eps_f * max(size_f, f_k) for LowIts iter, stop

% LowIts: No of iterations with low reduction before convergence

% SplineSmooth, SplineTol:
% Parameters which must be defined when computing a numerical approximation
% of the gradient and the Jacobian by use of the SPLINE Toolbox routines
% csaps.m and spaps.m:

% MajorIter: Maximal number of iterations for major problems 
% (Max Major Iterations)
% MinorIter: Maximal number of iterations for sub problems(Max Minor Iterations)

% MinorTolX: Subproblem convergence tolerance, in optimal sub-solution x
% Relative distance between successive x. //x_k+1 - x_k//. 

% PriFreq: Print frequency
% SummFreq: Summary frequency
% DiffInt: Difference Interval
% CentralDiff: Central difference Interval
% DerLevel: Derivative Level. Knowledge about gradients
% BigStep: Unbounded Step Size.
% BigObj: Unbounded Step Size.  Unbounded Objective.

% (7) LineAlg: Line search algorithm. 0 = quadratic, 1 = cubic, 3 = curved 
% (if available)

% (18) Initial step length. (Default 1 or less). 
% (21) LineSearch.sigma: Line search accuracy; 0 < sigma < 1
% sigma=0.9 inaccurate line search. sigma 0.1 accurate line search

% MODIFICATION LOG:
%
% 980918  hkh  Added line search parameters
% 980920  hkh  Added Penalty parameter for constrained problems
% 980922  hkh  Added bTol, eps_relf and LowIts. Also NOT_release_all.
% 981005  hkh  Added new field items, size_x,size_f,size_c,xTol,cTol,fTol,
%              bTol (new meaning), eps_dirg, eps_absf.
%              Deleted eps_fAbs,eps_fRel, eps_relf.
% 981010  hkh  Added new field subalg
% 981019  hkh  Added new field LineSearch.MaxIter, with default value 15
% 981026  hkh  Delete f_Low from optParam struct, put on Prob.f_Low
% 981027  hkh  Delete fields alg and subalg
% 981110  mbk  Add fields splineSmooth, default 0.4, and splineTol
%         hkh  Change eps2 to 100*eps
% 981210  hkh  Change to cubic linesearch as default
% 990910  hkh  Add parameter optParam.IterPrint
% 000709  hkh  Add parameter optParam.MinorIter
% 000910  hkh  Revision. Use SOL defaults, add some new fields.
% 010715  hkh  Added glbFast
% 010815  hkh  Added glcFast
% 011111  hkh  Added glcCluster and rbfSolve
% 020512  hkh  Adding bqpd
% 020622  hkh  Adding MIQPbb and changing parameters for bqpd as well
% 020630  hkh  Adding MINLPbb, filterSQP and PENSDP
% 020702  hkh  Avoid using strcmp, only use strcmpi
% 020921  hkh  Adding CPLEX
% 030116  hkh  Change name to Tlsqr, Add PENBMI
% 030126  hkh  Add PDCO, PDSCO. 
% 030126  hkh  Set WAIT and PRESOLVE before solvers, makes changes possible
% 030427  hkh  Add KS solver nlpql
% 040311  hkh  Avoid DiffInt and CentralDiff to be set as -900
% 040326  hkh  Set xTol lower for some solver, correct wrong comments about xTol
% 040528  hkh  Added XA, changed to 2000 for MaxIter for all LP,QP,MILP,MIQP
% 050308  hkh  Default splineSmooth -1 (let spaps determine optimal value)
% 050310  frhe Added glbDirect and glcDirect as aliases for glbFast and glcFast
% 060818  hkh  IMPORTANT CHANGE: eps_absf now small value close to eps
% 070221  hkh  Add arbfmip and arbf
% 080310  hkh  Change default iterations to 2000 from 500 for nlpSolve

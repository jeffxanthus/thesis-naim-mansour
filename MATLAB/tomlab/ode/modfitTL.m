% TOMLAB MODFIT Ordinary Differential Equations Solvers
%
% function Result = modfitTL(Prob)
%
% INPUT:
% Prob                    Problem structure in TOMLAB format.
%
% Prob.Solver.Name        Which MODFIT solver to use for optimization
%                         These are currently availible:
%                           dfnlp
%                           dn2gb   - default
%                           dlsmdf
%                           dfneld  - throws I/O recursion error in LNX
%
% Prob.SolverODE          Which MODFIT odesolver to use
%                         These are currently availible:
%                           dopri5  - default
%                           radau5
%                           ind-dir
%
% OUTPUT:
% Result     Structure with results (see ResultDef.m):

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Apr 21, 2005.  Last modified Jul 17, 2009.

function Result = modfitTL(Prob)

if nargin < 1
   error('modfitTL needs the Prob structure as input');
end

% All this should not be necessary if odeFit becomes ProbType
Prob.solvType = checkType('cls');
% Prob.PreSolve = 0;
% Prob.optParam.PreSolve = 0;
Prob = iniSolve(Prob,Prob.solvType,0,0);

Result=ResultDef(Prob);

Result.Solver='MODFIT';
Result.SolverAlgorithm=['solver ' Prob.Solver.Name ' and ODE-solver ' Prob.SolverODE '.'];

% Setup modfit with Tomlab Prob structure as input

% Temporary variables
E            = Prob.ODE.E;
nS           = Prob.ODE.nSeries;
Em           = length(E)/nS;
t            = Prob.ODE.tWant;
nt           = length(t);
rWht         = Prob.ODE.resWeight;

% Prepare measurements data
E = reshape(E,Em,nS);
if isempty(rWht)
   rWht = ones(Em,nS);
end
rWht = reshape(rWht,Em,nS);
NaNIdx   = find(isnan(E));
rWht(NaNIdx) = 0;

% Fitting criteria function
fit.npar    = Prob.N;
fit.x_0     = Prob.x_0;
fit.x_L     = Prob.x_L;
fit.x_U     = Prob.x_U;
fit.m_eval  = 'modfit_eval';

% Measurement data
meas.ntime  = nt;
meas.nconc  = 0;
meas.nmeas  = nS;
meas.time   = t(:);
meas.conc   = [];
meas.ymeas  = E;

if meas.nconc == 0
   meas.weight = rWht;
else
   meas.weight = repmat(rWht,[size(rWht) nconc]);
end

% Ordinary differential equations
ode.ntime   = nt;
ode.node    = length(Prob.ODE.Y0);
ode.ndae    = 0;

con         = [];
disc        = [];

par.params = ones(20,1)*NaN;

%       1    ioptp1 Number of iterations (default 120)
par.params(1)  = DefPar(Prob.optParam,'MaxIter',10000);

Prob.Solver.Name = DefPar(Prob.Solver, 'Name', 'dn2gb');
switch lower(Prob.Solver.Name)
  case 'dfnlp'
    par.method = 1;
%       2    ioptp2
% method = 1 maximum number of line search iterations (20)
%       4    opte1
% method = 1 final termination tolerance (1E-8)
%       5    opte2
    par.params(4) = DefPar(Prob.optParam, 'eps_f', NaN);
% Tolerance for chosen optimization algorithm:
% method = 1 expected size of residual
%       6    opte3
% Tolerance for chosen optimization algorithm:
%            not used
  case 'dn2gb'
    par.method = 3;
%       2    ioptp2
% method = 3 Maximum number of function calls
    par.params(2) = DefPar(Prob.optParam, 'MaxFunc', NaN);
%       4    opte1
% method = 3 relative accuracy in parameters
    par.params(4) = DefPar(Prob.optParam, 'xTol', NaN);
%       5    opte2
% Tolerance for chosen optimization algorithm:
% method = 3 relative error in function values
    par.params(5) = DefPar(Prob.optParam, 'fTol', NaN);
%       6    opte3
% Tolerance for chosen optimization algorithm:
%            not used
    case 'dslmdf'
    par.method = 4;
%       2    ioptp2
% method = 4 maximum number of function calls in search step
%       4    opte1
% method = 4 termination tolerance
    par.params(4) = DefPar(Prob.optParam, 'eps_f', NaN);
%       5    opte2
% Tolerance for chosen optimization algorithm:
% method = 4 initial steplength for search
%       6    opte3
% Tolerance for chosen optimization algorithm:
% method = 4 steplength reduction factor
  case 'dfneld'
    par.method = 5;
%       2    ioptp2
%            not used
%       4    opte1
% method = 5 accuracy in function values
    par.params(5) = DefPar(Prob.optParam, 'fTol', NaN);
%       5    opte2
% Tolerance for chosen optimization algorithm:
% method = 5 reflection factor
%       6    opte3
% Tolerance for chosen optimization algorithm:
% method = 5 contraction/expansion factor
    otherwise
    error('Illegal modfit-solver');
end

%       3    ioptp3 Not used
% ? par.params(4) = DefPar(Prob.optParam, 'eps_f' , 1e-8); % opte1
%       7    iodep1 Choice of odesolver
Prob.SolverODE = DefPar(Prob, 'SolverODE', 'dopri5');
switch lower(Prob.SolverODE)
  case 'dopri5'
    par.params(7) = 1;
%  case 'dop853'
%    par.params(7) = 2;
%  case 'odex'
%    par.params(7) = 3;
  case 'radau5'
    par.params(7) = 4;
%  case 'sdirk4'
%    par.params(7) = 5;
%  case 'seulex'
%    par.params(7) = 6;
  case 'ind-dir'
    par.params(7) = 11;
  case []
    par.params(7) = 1; % default
  otherwise
    error('Illegal modfit-odesolver');
end

%       8    iodep2 Only for DAE and implicit methods
%       9    iodep3 Approximate number of correct digits when gradients
%                   must be evaluated numerically using forward differences
%      10    iodep4 Jacobian bandwith for odesolvers RADAU5, SDIRK4, SEULEX
%      11    odee1  FINAL TERMINATION ACCURACY FOR ODE-SOLVER WITH RESPECT TO
%                   THE RELATIVE GLOBAL ERROR. (1E-7)
par.params(11)  = DefPar(Prob.ODE, 'relTol' , 1e-7);

%      12    odee2  FINAL TERMINATION ACCURACY FOR ODE-SOLVER WITH RESPECT TO
%                   THE ABSOLUTE GLOBAL ERROR. (1E-7)
par.params(12)  = DefPar(Prob.ODE, 'absTol' , 1e-7);

%      13    odee3  INITIAL STEPSIZE FOR SOLVING DIFFERENTIAL EQUATION
%                   (1E-4, only for DGEAR and RADAU)
par.params(13)  = DefPar(Prob.ODE, 'InitStep' , NaN);

Prob.NLPQL = DefPar(Prob, 'NLPQL', []);
%      14    nlpmi  Maximum number of iterations for NLPQL (50)
par.params(14)  = DefPar(Prob.NLPQL, 'maxit', NaN);
%      15    nlpac  Convergence tolerance for NLPQL (1E-12)
par.params(15)  = DefPar(Prob.NLPQL, 'acc', NaN);
%      16    nlpip  NLPQL print flag

%      17    norm   For method in [0,5], norm used for data fitting, where L1-
%                   or maxiumum norm are applicable only for simulation or DFNLP
%            norm=1   L1-NORM, I.E., SUM OF ABSOLUTE VALUES OF DEVIATIONS
%            norm=0/2 L2-NORM, I.E., SUM OF SQUARES OF DEVIATIONS
%            norm=3   LINF-NORM, I.E., MAXIMUM OF ABSOLUTE VALUES OF DEVIATIONS
%            norm   For method in [6,7], norm is design CRITERION SUBJECT
%                   TO INFORMATION MATRIX:
%            norm=1  -  SUM OF ABSOLUTE DIAGONAL VALUES MINIMIZED
%            norm=2  -  SMALLEST LARGEST EIGENVALUE MAXIMIZED
%            norm=3  -  CONDITION NUMBER MINIMIZED
%            norm=4  -  DETERMINANT MAXIMIZED
%            norm=5  -  AVERAGE VARIANCE MINIMIZED
%            norm=6  -  AVARAGE STANDARD DEVIATION MINIIMIZED, I.E.,
%                       AVERAGE LENGTH OF SIGNIFICANCE INTERVAL
par.params(17) = 2;     % Use L2-norm

%      18    numgra     NUMERICAL GRADIENT EVALUATION: (default 1)
%            numgra=-1  - ANALYTICAL DERIVATIVES PROVIDED BY THE USER
%            numgra=0/1 - FORWARD DIFFERENCES
%            numgra=2   - TWO-SIDED DIFFERENCES
%            numgra=3   - 5-POINT DIFFERENCE FORMULA
par.params(18) = 1;     % Use forward differences

%      19    nscale     SCALING (0-NO,1-STARTING VALUE(default),2-LOGARITHMIC)
par.params(19) = 0;     % Scale method = 1

%      20    isht       Shooting index (default 0)

prilev = 3;

outfile = ['Out.txt'];
scrfile = ['Scr.txt'];

[x, ifail, ymodel, g, f] = modfitmex(fit, meas, ode, con, disc, par,...
                                     Prob, prilev, outfile, scrfile);
Result.ODE.y = ymodel;
Result.Inform = ifail;
Result.x_0 = Prob.x_0;
Result.x_k = x;
Result.g_k = g;
Result.f_k = 0.5*(f(:)'*f(:));

Result=endSolve(Prob,Result);

% 050421  bkh  Written
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to
%              InitStep, tInit and tStop respectively
% 050429  bkh  Works for three solvers and three odesolvers
% 050705  med  Help updated
% 090717  med  Residual calculation updated

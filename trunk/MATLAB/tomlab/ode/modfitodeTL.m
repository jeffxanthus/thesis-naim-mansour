% TOMLAB MODFIT Ordinary Differential Equations Solver
%
% function Result = modfitTL(Prob, Solver)
%
% INPUT:
% Prob       Problem structure in TOMLAB format.
% Solver     Which MODFIT odesolver to use
%            These are currently availible:
%              dopri5   - default
%              radau5
%              ind-dir
%
% OUTPUT:
% Result     Structure with results (see ResultDef.m):

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Result = modfitodeTL(Prob, Solver)

if nargin < 2
   Solver = 'radau5';
   if nargin < 1 
      error('modfitodeTL needs the Prob structure as input');
   end
end

% Setup modfit for solving ode with Tomlab Prob structure as input

% Temporary variables
t              = Prob.ODE.tWant(:);
nt             = length(t);
Y0             = Prob.ODE.Y0;
nEq            = length(Y0);

% Necessary fields in Prob for modfit_eval
Prob.ODE.Y0Idx = DefPar(Prob.ODE, 'Y0Idx', []);
Prob.ODE.X = DefPar(Prob.ODE, 'X', 0);

% Fitting criteria function
fit         = [];
fit.npar    = DefPar(Prob, 'N'  , length(Prob.ODE.X));
fit.x_0     = Prob.ODE.X;
fit.x_L     = -inf*ones(fit.npar,1);
fit.x_U     = inf*ones(fit.npar,1);
fit.m_eval  = 'modfit_eval';

% Measurement data
meas = [];
meas.ntime  = nt;
meas.nconc  = 0;
meas.time   = t(:);
meas.conc   = [];

% Ordinary differential equations
ode.ntime   = nt;
ode.node    = nEq;
ode.ndae    = 0;

con         = [];
disc        = [];

par.params = ones(20,1)*NaN;

% No optimization, only evalutation
par.method = 0;

% par.params(i) for i = 1:20

%       i
%       1    ioptp1 Confidence level (1/5/10)
%       2    ioptp2 N/A
%       3    ioptp3 Not used
%       4    opte1  Round-off tolerance for determining rank of cov matrix
%       5    opte2  N/A
%       6    opte3  N/A
%       7    iodep1 Choice of odesolver
switch lower(Solver)
  case 'dopri5'
    par.params(7) = 1;
%  case 'dop853'
%    par.params(7) = 2;
%  case 'odex'
%    par.params(7) = 3;
  case 'radau5'        % implicit methods not availible yet
    par.params(7) = 4;
%  case 'sdirk4'
%    par.params(7) = 5;
%  case 'seulex'
%    par.params(7) = 6;
  case 'ind-dir'
    par.params(7) = 11;
  case {'radau5', 'sdirk4', 'seulex'}
    fprintf('Implicit odesolvers not availible yet');
    error('Illegal choice of modfit-odesolver');
  case []
    par.params(7) = 1; % default
  otherwise      
    error('Illegal choice of modfit-odesolver');
end
%par.params(7) = NaN;
%       8    iodep2 Only for DAE and implicit methods
%       9    iodep3 Approximate number of correct digits when gradients
%                   must be evaluated numerically using forward differences
%      10    iodep4 Jacobian bandwith for odesolvers RADAU5, SDIRK4, SEULEX
%      11    odee1  FINAL TERMINATION ACCURACY FOR ODE-SOLVER WITH RESPECT TO
%                   THE RELATIVE GLOBAL ERROR. (1E-7)
par.params(11)  = DefPar(Prob.ODE, 'relTol' , 1e-12);

%      12    odee2  FINAL TERMINATION ACCURACY FOR ODE-SOLVER WITH RESPECT TO 
%                   THE ABSOLUTE GLOBAL ERROR. (1E-7)
par.params(12)  = DefPar(Prob.ODE, 'absTol' , 1e-12);

%      13    odee3  INITIAL STEPSIZE FOR SOLVING DIFFERENTIAL EQUATION 
%                   (1E-4, only for DGEAR and RADAU)
%par.params(13)  = DefPar(Prob.ODE, 'InitStep' , NaN);

%Prob.NLPQL = DefPar(Prob, 'NLPQL', []);
%      14    nlpmi  Maximum number of iterations for NLPQL (50)
%par.params(14)  = DefPar(Prob.NLPQL, 'maxit', NaN);
%      15    nlpac  Convergence tolerance for NLPQL (1E-12)
%par.params(15)  = DefPar(Prob.NLPQL, 'acc', NaN);

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

for i = 5:5:ceil(nEq/5)*5
   if i < nEq
      iEq = i;
   else
      iEq = nEq;
   end
   Prob.ODE.Eeq = (i-4:iEq)';
   meas.nmeas  = iEq-i+5;
   meas.ymeas  = NaN*ones(nt,meas.nmeas);
   meas.weight = ones(nt,meas.nmeas);
   [x, ifail, ymodel(:,Prob.ODE.Eeq), g, f] = modfitmex(fit, meas, ode, con, disc, par, ...
                    Prob, prilev, outfile, scrfile);                
end

Result.ODE.y = ymodel;
Result.Inform = ifail;

% MODIFICATION LOG:
%
% 050421  bkh  Written
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050428  bkh  Works for three methods
% 050705  med  Help updated
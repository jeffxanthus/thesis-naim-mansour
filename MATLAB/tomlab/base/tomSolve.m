% function Result = tomSolve(Solver,Prob)
%
% tomSolve solves a sub problem without any check on the Prob structure
%
% It is intended for use by other solvers when solving a subproblem like a
% QP, LP, Dual LP or FP (feasible point) problem
%
% Global variables are saved before the solver call, and restored afterwards
%
% It could also be used for recursive solutions of user problems
% or control loops when speed is a demand.
%
% INPUT PARAMETERS
% Solver    Name of the solver
% Prob      Input structure, feeded to the solver
%
% OUTPUT PARAMETERS
% Result    Output result structure, feeded back to the caller.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Nov 6, 2000.    Last modified May 8, 2009.

function Result = tomSolve(Solver,Prob)

if nargin < 2
   error('SolveQP needs input structure Prob');
end

global GlobalLevel

if isempty(GlobalLevel)
   GlobalLevel=1;
else
   GlobalLevel=GlobalLevel+1;
end
Level=GlobalLevel;
globalSave(Level);

switch lower(Solver)
   case {'minos'}
      Result = minosTL(Prob);
   case {'milpsolve'}
      Result = milpSolveTL(Prob);
   case {'npsol'}
      Result = npsolTL(Prob);
   case {'nlssol'}
      Result = nlssolTL(Prob);
   case {'snopt','snopt6','snopt7'}
      Result = snoptTL(Prob);
   case {'snopt8'}
      Result = snopt8TL(Prob);
   case {'ucsolve'}
      Result = ucSolve(Prob);
   case {'consolve'}
      Result = conSolve(Prob);
   case {'nlpsolve'}
      Result = nlpSolve(Prob);
   case {'clssolve'}
      Result = clsSolve(Prob);
   case {'glbsolve'}
      Result = glbSolve(Prob);
   case {'glbfast'}
      Result = glbFast(Prob);
   case {'glcsolve'}
      Result = glcSolve(Prob);
   case {'glcfast'}
      Result = glcFast(Prob);
   case {'glccluster'}
      Result = glcCluster(Prob);
   case {'rbfsolve'}
      Result = rbfSolve(Prob);
   case {'arbfmip'}
      Result = arbfmip(Prob);
   case {'mipsolve'}
      Result = mipSolve(Prob);
   case {'multimin'}
      Result = multiMin(Prob);
   case {'strustr'}
      Result = sTrustr(Prob);
   case {'filtersqp'}
      Result = filterSQPTL(Prob);
   case {'nlpql','nlpqlp'}
      Result = nlpqlTL(Prob);
   case {'pdco'}
      Result = pdcoTL(Prob);
   case {'pdsco'}
      Result = pdscoTL(Prob);
   case {'oqnlp'}
      Result = oqnlpTL(Prob);
   case {'msnlp'}
      Result = msnlpTL(Prob);   
   case {'knitro'}
      Result = knitroTL(Prob);
   case {'conopt'}
      Result = conoptTL(Prob);
   case {'minlpbb'}
      Result = minlpBBTL(Prob);
   case {'fmincon'}
      Result = opt20Run('fmincon',Prob);
   case {'fminsearch'}
      Result = opt20Run('fminsearch',Prob);
   case {'fminunc'}
      Result = opt20Run('fminunc',Prob);
   case {'lsqcurvefit'}
      Result = opt20Run('lsqcurvefit',Prob);
   case {'lsqnonlin'}
      Result = opt20Run('lsqnonlin',Prob);
   case {'lsqnonneg'}
      Result = opt20Run('lsqnonneg',Prob);
   case {'constr'}
      Result = opt15Run('constr',Prob);
   case {'fmins'}
      Result = opt15Run('fmins',Prob);
   case {'fminu'}
      Result = opt15Run('fminu',Prob);
   case {'leastsq'}
      Result = opt15Run('leastsq',Prob);
   case {'ego'}
      Result = ego(Prob);
   case {'ego05'}
      Result = ego05(Prob);
   case {'quadprog'}
      Result = opt20Run('quadprog',Prob);
   case {'linprog'}
      Result = opt20Run('linprog',Prob);
   case {'miqpbb'}
      Result = miqpBBTL(Prob);
   case {'cutplane'}
      Result = cutplane(Prob);
   case {'xpress-mp'}
      Result = xpressTL(Prob);
   case {'cplex'}
      Result = cplexTL(Prob);
   case {'gurobi'}
      Result = gurobiTL(Prob);
   case {'xa'}
      Result = xaTL(Prob);
   case {'bqpd'}
      Result = bqpdTL(Prob);
   case {'pensdp'}
      Result = pensdpTL(Prob);
   case {'penbmi'}
      Result = penbmiTL(Prob);
   case {'qpsolve'}
      Result = qpSolve(Prob);
   case {'lpsimplex'}
      Result = lpSimplex(Prob);
   case {'dualsolve'}
      Result = DualSolve(Prob);
   case {'lsei'}
      Result = lseiTL(Prob);
   case {'tlsqr'}
      Result = TlsqrTL(Prob);
   case {'lssol'}
      Result = lssolTL(Prob);
   case {'sqopt'}
      Result = sqoptTL(Prob);
   case {'gp', 'coplgp'}
      Result = coplgpTL(Prob);
   case {'qpopt'}
      Result = qpoptTL(Prob);
   case {'lpopt'}
      Result = lpoptTL(Prob);
   case {'lp-minos'}
      Result = minoslpTL(Prob);
   case {'qld'}
      Result = qldTL(Prob);
   case {'qp-minos'}
      Result = minosqpTL(Prob);
   case {'qp'}
      Result = opt15Run('qp',Prob);
   case {'lp'}
      Result = opt15Run('lp',Prob);
   otherwise
      Result = tomRun(Solver,Prob);
end

globalGet(Level);
GlobalLevel=GlobalLevel-1;

% MODIFICATION LOG:
%
% 001106  hkh  Written
% 010717  hkh  Added glbFast and Xpress-MP
% 010815  hkh  Added glcFast
% 011104  hkh  Added glcCluster
% 020526  hkh  Adding Dundee QP solver bqpd
% 020621  hkh  Adding Dundee MIQP solver MIQPbb
% 020630  hkh  Adding Dundee solvers MINLPbb, filterSQP; and PENSDP
% 020702  hkh  Adding CPLEX
% 030116  hkh  Change to Tlsqr, add PENBMI
% 030123  hkh  Add pdco, pdsco
% 030129  hkh  Remove empty varargin in call to opt15Run and opt20Run
% 040101  hkh  Add nlpql
% 040103  hkh  Direct call to all TL files
% 040312  ango Spelling of qpoptTL, lpoptTL fixed
% 041216  med  OQNLP and MSNLP added
% 050112  med  Added LSGRG2, KNITRO and CONOPT, XA
% 050214  frhe Added BOEING solvers
% 050323  med  Added LPSOLVE and moved lpSolve to lpSolve2
% 050330  ango lpSolve2 is now lpSimplex
% 050423  ango lpSolve changed to milpSolve
% 050602  hkh  Add coplgp (also possible to call as gp)
% 050605  hkh  Add snopt7, snopt8, and use snopt and snopt6 to call snopt6
% 060201  ango Add ego05
% 060212  hkh  Add arbfmip
% 060715  hkh  Add multiMin
% 080607  hkh  Change order, last LP, QP, MILP etc.
% 081030  ango snopt7 replaces snopt,snopt6
% 090508  med  gurobi added
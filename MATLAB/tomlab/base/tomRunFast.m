% tomRunFast - Sub problem driver routine for TOMLAB
%
% call with:
%
%   function Result = tomRunFast(Solver, Prob);
%
% tomRun assumes first argument is a string
% tomRun assumes second argument is a structure
%
% Solver     The name of the solver that should be used to optimize the problem.
%
% Prob       Problem structure.
%
% Result Structure with optimization results

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jun 7, 2008.    Last modified Sept 11, 2009.

function Result = tomRunFast(Solver, Prob)

global n_f n_g n_H    % Count of function, gradient, Hessian evaluations
global n_c            % Count of constraint, constr.grad, and 2nd der evals
global n_r n_J        % Count of residual, Jacobian and 2nd der evaluations
global GlobalLevel optType probType
GlobalLevel = [];     % Initialize to zero depth for recursion
probType    = [];
optType     = [];
NLLSVars = 0;
GlobVars = 0;

switch lower(Solver)
 case 'nlpsolve'
   Result = nlpSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;
 case 'consolve'
   Result = conSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;
 case 'strustr'
   Result = sTrustr(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
 case 'clssolve'
   Result = clsSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;
 case 'glbsolve'
   Result= glbSolve(Prob);
 case 'glbfast'
   Result= glbFast(Prob);
 case 'glbdirect'
   Result= glbDirectTL(Prob);
 case 'ego'
   Result= ego(Prob);
 case 'ego05'
   Result= ego05(Prob);
 case 'glcsolve'
   Result= glcSolve(Prob);
 case 'glcfast'
   Result= glcFast(Prob);
 case 'glcdirect'
  Result= glcDirectTL(Prob);
 case 'glccluster'
   Result= glcCluster(Prob);
 case 'minlpsolve'
   Result = minlpSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   NLLSVars = 1;
   GlobVars = 1;
 case 'rbfsolve'
   Result= rbfSolve(Prob);
 case 'arbfmip'
   Result= arbfmip(Prob);
 case 'arbf'
   Result= arbf(Prob);
 case 'dynrbf'
   Result= dynrbf(Prob);
 case 'milpsolve'
   Result = milpSolveTL(Prob);
 case 'multimin'
   Result = multiMin(Prob);
 case 'multiminlp'
   Result = multiMINLP(Prob);
 case 'sepminlp'
  Result = sepMINLP(Prob);
 case 'lpsimplex'
   Result = lpSimplex(Prob);
 case 'qpsolve'
   Result = qpSolve(Prob);
 case 'mipsolve'
   Result = mipSolve(Prob);
 case 'cutplane'
   Result = cutplane(Prob);
 case 'ucsolve'
   Result = ucSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;
 case 'lpopt'
   Result = lpoptTL(Prob);
 case 'qpopt'
   Result = qpoptTL(Prob);
 case {'sqopt','sqopt7'}
   Result = sqoptTL(Prob);
 case 'lssol'
   Result = lssolTL(Prob);
 case 'nlssol'
   Result = nlssolTL(Prob);
 case 'minos'
   Result = minosTL(Prob);
   NLLSVars = 1;
 case 'lp-minos'
   Result = minoslpTL(Prob);
 case 'qp-minos'
   Result = minosqpTL(Prob);
 case 'npsol'
   Result = npsolTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case {'snopt','snopt6','snopt7'}
   Result = snoptTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'snopt8'
   Result = snopt8TL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'tlsqr'
   Result = TlsqrTL(Prob);
 case 'qld'
   Result = qldTL(Prob);
 case 'lsei'
   Result = lseiTL(Prob);
 case 'bqpd'
   Result = bqpdTL(Prob);
 case 'miqpbb'
   Result = miqpBBTL(Prob);
 case {'filtersqp', 'filsqp'}
   Result = filterSQPTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'minlpbb'
   Result = minlpBBTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'pensdp'
   Result = pensdpTL(Prob);
 case 'penbmi'
   Result = penbmiTL(Prob);
 case {'nlpql', 'nlpqlp'}
   Result = nlpqlTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case {'dfnlp', 'dfnlpd'}
   Result = dfnlpTL(Prob);
 case {'nlpjob'}
   Result = nlpjobTL(Prob);
 case {'misqp'}
   Result = misqpTL(Prob);
 case 'pdco'
   Result = pdcoTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'pdsco'
   Result = pdscoTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'oqnlp'
   Result  = oqnlpTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'msnlp'
   Result  = msnlpTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;   
 case 'lsgrg2'
      Result  = lsgrg2TL(Prob);
   NLLSVars = 1;
   GlobVars = 1;   
 case {'knitro'}
      Result  = knitroTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;   
 case {'conopt'}
      Result = conoptTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case {'xpress-mp','xpress','xpressmp'}
      Result = xpressTL(Prob);
 case 'cplex'
      Result = cplexTL(Prob);
 case 'cplex11'
      Result = cplex11TL(Prob);
 case 'gurobi'
      Result = gurobiTL(Prob);
 case 'xa'
      Result = xaTL(Prob);
 case 'constr'
   Result=opt15Run('CONSTR',Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'fmins'
   Result=opt15Run('FMINS',Prob);
 case 'leastsq'
   Result=opt15Run('LEASTSQ',Prob);
 case 'lp'
   Result=opt15Run('LP',Prob);
 case 'qp'
   Result=opt15Run('QP',Prob);
 case 'fminu'
   Result=opt15Run('FMINU',Prob);
 case 'fminunc'
   Result=opt20Run('FMINUNC',Prob);
 case 'fmincon'
   Result=opt20Run('FMINCON',Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'fminsearch'
   Result=opt20Run('FMINSEARCH',Prob);
 case 'lsqnonlin'
   Result=opt20Run('LSQNONLIN',Prob);
 case 'lsqlin'
   Result=opt20Run('LSQLIN',Prob);
 case 'linprog'
   Result=opt20Run('LINPROG',Prob);
 case 'quadprog'
   Result=opt20Run('QUADPROG',Prob);
 case {'gp', 'coplgp'}
   Result=coplgpTL(Prob);
 case 'geno'
   Result=genoTL(Prob); 
 case 'lgo'
   Result=lgoTL(Prob);
 case 'slssolve'
   Result=slsSolve(Prob,0);
 case 'infsolve'
   Result=infSolve(Prob,0);
 case 'inflinsolve'
   Result=infLinSolve(Prob,0);
 case 'goalsolve'
   Result=goalSolve(Prob,0);
 case 'l1solve'
   Result=L1Solve(Prob,0);
 case 'l1linsolve'
   Result=L1LinSolve(Prob,0);
 case 'linratsolve'
   Result=linRatSolve(Prob,0);
 case 'expsolve'
   Result=expSolve(Prob,0);
 otherwise
   fprintf('Solver %s',Solver);
   fprintf(' NOT found\n');
   error('Illegal solver algorithm!')
end

if isempty(Result), EmptyResult(Solver), end

if any(Prob.probType==[4 5 6])
   if ~isempty(Result)
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end

% HKH CHECK THIS CAREFULLY, only used in mex routines before
if NLLSVars
   if any(Prob.probType==[4 6 11])
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end
if GlobVars
   if any(Prob.probType==[1 2 3 10])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   elseif any(Prob.probType==[1 9])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
   elseif any(Prob.probType==[4 6 11])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   end
end

function EmptyResult(Code)
fprintf('\n\n')
fprintf('tomRun: ');
fprintf('ERROR! %s solver returned empty Result',Code);
fprintf('\n\n')
fprintf('Run startup to set paths\n');
fprintf('\n\n')

% MODIFICATION LOG:
%
% 080607  hkh  Written, simplified version of tomRun
% 081112  ango sqopt7 replaces sqopt
% 090508  med  gurobi added
% 090911  hkh  multiMINLP added

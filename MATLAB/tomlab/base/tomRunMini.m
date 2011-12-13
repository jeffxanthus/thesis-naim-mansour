% tomRunMini - Silent driver routine for TOMLAB recursive calls
%
% Call with:
%
%   function Result = tomRun(Solver, Prob)
% 
% License check is already assumed to be done

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Jun 7, 2008.    Last modified May 8, 2009.

function Result = tomRunMini(Solver,Prob)

switch lower(deblank(Solver))
 case 'milpsolve'
   Result = milpSolveTL(Prob);
 case 'lpsimplex'
   Result = lpSimplex(Prob);
 case 'qpsolve'
   Result = qpSolve(Prob);
 case 'lpopt'
   Result = lpoptTL(Prob);
 case 'qpopt'
   Result = qpoptTL(Prob);
 case {'sqopt','sqopt7'}
   Result = sqoptTL(Prob);
 case 'lssol'
   Result = lssolTL(Prob);
 case 'lp-minos'
   Result = minoslpTL(Prob);
 case 'qp-minos'
   Result = minosqpTL(Prob);
 case 'tlsqr'
   Result = TlsqrTL(Prob);
 case 'qld'
   Result = qldTL(Prob);
 case 'lsei'
   Result = lseiTL(Prob);
 case 'bqpd'
   Result = bqpdTL(Prob);
 case {'barqp'}
   Result = barqpTL(Prob);
 case {'sprqp'}
   Result = sprqpTL(Prob);
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
 otherwise
   error('Illegal subsolver used');
end

% MODIFICATION LOG:
%
% 080607  med  Written
% 090508  med  gurobi added
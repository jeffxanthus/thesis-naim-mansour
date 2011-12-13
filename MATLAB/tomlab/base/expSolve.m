% function Result = expSolve(Prob, PriLev)
%
% Solve exponential fitting problems for given number of terms p
%
% expSolve calls Result = tomRun(SolverL2,Prob,PriLev);
%
% Also does statistical analysis of the solution calling StatLS.
% If statistical information is not required, tonRun can be called
% directly.
%
% INPUT:
% Solver  Name of solver to use. If empty, TOMLAB selects dependent on license
% Prob    Problem created with expAssign
% PriLev  Print level in tomRun call
%
% Prob.SolverL2  The subsolver to use
%
% OUTPUT:
% Result  Structure. Special fields in Result:
%
% Result.LS            Result of statistical analysis, see help StatLS
% Result.ExpFit.lambda Optimal intensities
% Result.ExpFit.alpha  Optimal weights
% Result.ExpFit.beta   [] unless eType == 4
%
% Note if SepAlg == 1:
%      The result from the optimization will only display the intensities
%      lambda, and the problem size Prob.N (in Result.Prob.N as output) is p.
%      Optimal alpha (and beta) values are found in the Result.ExpFit
%      structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Aug 11, 2002.    Last modified Feb 2, 2006.

function Result = expSolve(Prob, PriLev)

if nargin < 2
   PriLev = 0;
   if nargin < 1
      error('Input Prob must be given');
   end
end

global wLS LS_A

if isfield(Prob, 'SolverL2')
   Solver = Prob.SolverL2;
else
   Solver = [];
end

if isempty(Solver), Solver  = GetSolver('exp',0,0); end

if strcmpi('nlssol',Solver)
   % Special parameters for least squares problems are:
   % JTJ Initial Hessian (often best to have as true)
   if Prob.SOL.optPar(47) < -10
      Prob.SOL.optPar(47) = 1;  % Default 1, other unit Hessian
   end

   % RESET Frequency . When Gauss-Newton works fine, often for small
   % residual problems, then one may raise this value
   if Prob.SOL.optPar(48) < -10
      Prob.SOL.optPar(48) = 2;  % Default 2, Reset each 2nd step
   end
end

Result  = tomRun(Solver,Prob,PriLev);

% Parameter statistics, get Jacobian from Result structure
p = Prob.ExpFit.p;
eType = Prob.ExpFit.eType;

if Prob.LS.SepAlg
   % Determine alpha by NNLS. Jz is weighted by wLS
   % LS_A is E matrix of exponential expression
   [alpha,beta,Jz,wLS]=expLS(Result.x_k,LS_A,Prob);
   % Will result in wrong Jacobian for full dimension, because of projection
   Result.ExpFit.lambda = Result.x_k;
   Result.ExpFit.alpha  = alpha;
   Result.ExpFit.beta   = beta;
   Prob.LS.SepAlg = 0;
   x_k = [Result.x_k;alpha;beta];
   Prob.N = length(x_k);
   r_k = nlp_r(x_k,Prob);
   J_k = nlp_J(x_k,Prob);
   Result.LS = StatLS(x_k, r_k, J_k);
else
   Result.ExpFit.lambda=Result.x_k(1:p);
   Result.ExpFit.alpha=Result.x_k(p+1:2*p);
   if eType == 4
      Result.ExpFit.beta=Result.x_k(2*p+1:3*p);
   else
      Result.ExpFit.beta=[];
   end
   Result.LS = StatLS(Result.x_k, Result.r_k, Result.J_k);
end

% MODIFICATION LOG:
%
% 041123  hkh  Change call to tomRun
% 050802  hkh  Add eType as input
% 050802  hkh  Handle SepAlg correctly
% 060131  med  Assign code moved to expAssign
% 060201  med  Call changed to standard format, SolverL2 added
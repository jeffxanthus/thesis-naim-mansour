% TOMLAB GP Geometric Programming Solver
%
% function Result = coplgpTL(Prob)
%
% COPLGP solves geometric programming problems of the type:
%
%    min g0(t)
%
%    subject to
%
%    gk(t) <= 1, k =1,2...p
%
%    t(i)>= 0
%
% where
%
%    g0(t) = sum{j=1,n(0)} [c(j)*t(1)^a(j,1)...t(m)^a(j,m)]
%
%    gk(t) = sum{j=n(k-1)+1,n(k)} [c(j)*t(1)^a(j,1)...t(m)^a(j,m)]
%    for k = 1,2,...,p
%
%
% The COPLGP program solves posynomial GP problems. The dual solved is:
%
%        max    prod( (c(j)/x(j))^x(j)) * prod(lambda(k)^lambda(k))
%         x
%        s/t    sum(x(i)) = 1,       for i=1,2,...,n(0)
%               sum(x(j)*a(j,i)) = 0 for i=1,2,...,m
%               x(j) >= 0, j=1,2,...n(p)
%
% where lambda(k) = sum(x(j),j=n(k-1)+1:n(k)), k=1,2,...p
%
%
% The problem is formulated on the dual form, see gpAssign and the TOMLAB
% /GP manual for more information and examples.
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%  GP.A         Linear constraint matrix for dual problem.
%  GP.coef      Posynomial coefficient vector.
%  GP.nterm     Number of terms in objective and constraints.
%  GP.moremem   Extra memory to allocate. Default 300*n + 300*m,
%               where n is the total number of terms and m is the
%               number of variables. If this is not enough the
%               interface will double the memory until it is enough
%               and try to solve it. The final value of moremem
%               will be returned in Result.GP.lackmem. If you solve
%               many similar problems in a sequence, set
%               Prob.GP.moremem to Result.GP.lackmem in order to
%               avoid many unnecessary memory allocations.
%
%  PriLevOpt    Print level in the solver. See GP.options.PRILEV.
%
%  GP           Structure with GP solver specific fields.
%
%    PrintFile  Name of file to print progress information and
%               results to. The output written to file is always
%               with a print level of 2.
%               Default: Empty, none that is.
%
%  GP.options   Structure with special fields for the GP solver:
%
%      PRILEV   Print level of GP solver. Has precedence over
%               Prob.PriLevOpt.
%                 ==0   silent
%                 >=1   normal output with iteration log
%                 >=2   more verbose output
%                 >=10  memory debug output
%                 >=100 extreme debug output
%               Default: 0
%      ITLM     Iteration limit.
%               Default: 100
%      TOLX     Zero tolerance on the dual variables.
%               Default: 1e-14
%      TOLZ     Zero tolerance on the constraint values.
%               Default: 1e-14
%      EPSP     Feasibility tolerance on primal problem.
%               Default: 1e-6
%      EPSD     Feasibility tolerance on dual problem.
%               Default: 1e-6
%      EPSC     Feasibility tolerance on the linear dual constraints.
%               Default: 1e-9
%      TOLPIV   Pivoting zero tolerance.
%               Default: 1e-14
%      FRAC     Fractional decrease of steplength.
%               Default: 0.9995
%      BIG      Value to treat as infinite.
%               Default: 1e20
%
% -----------------------------------------------------------------------------
%
% OUTPUT:
% Result        Structure with optimization results.
%
%   f_k         Function value at optimum.
%
%   x_k         Solution vector.
%
%   xState      State of variables. Free == 0; On lower == 1; On upper == 2;
%               Fixed == 3;
%
%   ExitFlag    Exit status.
%   ExitText    Exit text from solver.
%   Inform      GP information parameter. See Result.ExitText
%
%   Solver           Name of the solver (GP).
%   SolverAlgorithm  Description of the solver.
%
%   GP.lackmem  See description for input Prob.GP.moremem.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written May 19, 2004.   Last modified Jun 8, 2008.

function Result = coplgpTL(Prob)

Prob.solvType = checkType('gp');

Prob = iniSolveMini(Prob);
Result = ResultDef(Prob);
Result.Solver = 'GP';
Result.SolverAlgorithm='COPLGP - Geometric Programming Solver';

GP             = DefPar(Prob,'GP',[]);
moremem        = DefPar(GP,'moremem',0);
options        = DefPar(GP,'options',[]);
PrintFile      = DefPar(GP,'PrintFile','');
PriLev         = DefPar(options, 'PRILEV', DefPar(Prob, 'PriLevOpt', ...
                                                  0));
options.PRILEV = PriLev;

n              = size(Prob.GP.A, 2);

[Inform, x_k, f_k, lackmem] = coplgp(Prob.GP.A,Prob.GP.coef,Prob.GP.nterm,options,PrintFile,moremem);

Result.GP.lackmem  = lackmem;

Result.xState = zeros(n, 1);
Result.xState(x_k <= 1e-8) = 1;

Result.x_k = x_k;
Result.f_k = f_k;

switch(Inform)
 case 0,    ExitText = 'Optimal solution obtained'; ExitFlag = 0;
 case 1,    ExitText = 'Primal is unbounded'; ExitFlag = 2;
 case 2,    ExitText = 'Dual is unbounded'; ExitFlag = 2;
 case 3,    ExitText = 'Too many iterations'; ExitFlag = 1;
 case 4,    ExitText = 'Numerical difficulties'; ExitFlag = 3;
 case 5,    ExitText = 'Infeasible detected'; ExitFlag = 4;
 case 6,    ExitText = 'Insufficient space'; ExitFlag = 6;
 case 7,    ExitText = 'Error in input data'; ExitFlag = 10;
 otherwise, ExitText = 'Unknown solver status'; ExitFlag = 999;
end

Result.Inform = Inform;
Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG
%
% 050519 med  File created
% 050526 frhe Documentation modified
% 050526 frhe Code adapted to work with mex
% 050527 frhe Updated problem formulation and help text
% 050528 med  xState defined with a tolerance of 1e-8
% 050601 med  solvType set dynamically
% 050602 hkh  Change name to coplgp and coplgpTL
% 080606 med  Switched to iniSolveMini
% 08060i hkh  Switched to endSolveMini

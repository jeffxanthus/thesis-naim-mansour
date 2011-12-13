% TOMLAB miqpBB MIQP/MILP, LP/QP Solver
%
% function [Inform,x_k,Obj,Iter] =
%         miqpbb(A, bl, bu, IntVars, Priority, Func,...
%                mlp, kmax, stackmax, optPar,...
%                PriLev, PrintFile, Prob, moremem)
%
%
% miqpBB solves the following mixed integer (linear or quadratic)
% programming problem (MILP, MIQP, LP, QP):
%
%   minimize   0.5 * x'*F*x + c'x
%      x
%
%  subject to     bl <= (  x ) <= bu
%                       ( Ax )
%
%  In addition, some or all x may be integer valued.
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints).
%   bl,bu has dimension n+m, with simple bounds in the n first elements.
%
%   c has dimension n.
%   F is a n x n symmetric matrix, sparse or dense.
%   If F is empty, an LP or MILP problem is solved.
%
% ------------------------------------------------------------------
%
% INPUT:
%
%
% [c A']    Linear constraint matrix, dense or sparse n x (m+1)
%           matrix. miqpBB requires the transpose of the constraint
%           matrix, and also with the linear part of the objective
%           function as the first column.
%
% bl,bu     Lower and upper bounds on variables and constraints. Length
%           must be n+m where the first n elements are simple bounds
%           on the variables.
%
% IntVars   Vector with integer variable indices.
%
% Priority  Priorities for the integer variables. Length must the
%           same as that of IntVars.
%
% Func      Name of MATLAB callback function that performs the
%           Hessian - vector multiplication F*x. A standard routine
%           is supplied in tomlab/lib/HxFunc.m, using the Prob.QP.F
%           matrix. If the user for some reason wants to write his
%           own callback function, it must take arguments similar
%           to those of HxFunc.m. The second argument nState is
%           always 0.0 in the current version of the solver.
%
% mlp       Maximum level parameter for resolving degeneracy in
%           BQPD which is used as sub-problem solver. If empty, the
%           MEX interface sets mlp to m, the number of constraints.
%
% kmax      Maximum dimension of reduced space. Default (and
%           maximum) value is n, the number of variables.
%
% stackmax  Size of the stack storing information during the
%           tree-search. Default value if empty: 5000.
%
%
% optPar    Vector of optimization parameters. If -999, set to default
%           Length from 0 to 20 allowed. The following elements are
%           used by miqpBB:
%
% optPar(1)  iprint 0    Print level in miqpbb
% optPar(2)  tol   1E-10 Relative accuracy in solutio
%            (Prob.optParam.eps_x)
% optPar(3)  emin  1.0   1.0 Use cscale (constraint scaling) 0.0 no scaling
% optPar(4)  sgnf  5E-4  Max rel error in two numbers equal in exact arithmetic
% optPar(5)  nrep  2     Max number of refinement steps
% optPar(6)  npiv  3     No repeat if no more than npiv steps were taken
% optPar(7)  nres  2     Max number of restarts if unsuccessful
% optPar(8)  nfreq 500   The max interval between refactorizations
%
%
% optPar(12) epsilon     Tolerance used for x value tests. DEF 1E-5
%            (Prob.optParam.eps_f)
% optPar(13) MIopttol    Tolerance used for function value tests 1E-4
%
% optPar(14) fIP. Upper bound on f(x). Only consider solutions < fIP - epsilon
%            Prob.MIP.fIP. DEFAULT infty = optPar(9) = 1E20
%
% optPar(15) timing      = 1, Use timing, = 0 no timing (default)
% optPar(16) max_time    Maximal time allowed for the run in seconds
%            Default 4E3, i.e. 66 minutes
%
% optPar(17) branchType  Branch on variable with highest priority. If tie:
%            = 1. Variable with largest fractional part, among those branch
%            on the variable giving the largest increase in the objective
%            = 2. Tactical Fletcher (PMO) branching. The var that solves
%            max(min(e+,e-)) is chosen. The problem than corresponding to
%            min(e+,e-) is placed on the stack first.
%            = 3. Tactical branching, Padberg/Rinaldi,91, Barahona et al.,89
%            (i) Choose the branching variable the one that most violates the
%            integrality restrictions. i.e. find  max(i){min(pi+,pi-)}
%            pi+ = int(x(i)+1) - x(i) , pi- = x(i) - int(x(i))
%            (ii) among those branch on the variable that gives the greatest
%            increase in the obj. function (iii) Finally a LOWER BOUND is
%            computed on the branched problems using the bounding method of
%            Fletcher and Leyffer (dual active set step) DEFAULT = 1
%
% optPar(18) ifsFirst If 1, then only search for first ifs (ifail=6),DEFAULT 0
%
% optPar(19) infty       Real value for infinity  (default 1E20)
%
%
% PriLev     Print level in the MEX interface:
%            0 = off, 1 = only result is printed, 2 = result and
%            intermediate steps are printed. scalar information, 3 = verbose).
%
% PrintFile  Name of print file. Amount/print type determined by
%            PriLev parameter. Default name miqpbbout.txt.
%
% Prob       The Tomlab problem description structure. This is a
%            necessary argument if the standard HxFunc.m callback
%            routine is used. HxFunc uses Prob.QP.F to calculate
%            the Hessian*vector multiplication.
%
% moremem    Scalar or 2x1-vector giving values for extra work
%            memory allocation. If scalar, the value given is added
%            to both the INTEGER and REAL workspaces.
%            If a vector is given, the first element controls the REAL
%            workspace increase and the second the INTEGER
%            workspace. Set one or both elements to values <0 for
%            problem dependent memory increases.
%
% -----------------------------------------------------------------
%
% OUTPUT
%
% ifail      Status code: the following values are defined:
%            0 - Solution found
%            1 - Error in parameters for BQPD
%            2 - Unbounded QP encountered
%            3 - Stack overflow - no integer solution found
%            4 - Stack overflow - some integer solution found
%            5 - Integer infeasible
%            7 - Infeasible root problem
%
% x_k        The solution vector, if any found. If ifail is other
%            than 0 or 4, the contents of x is undefined.
%
% Obj        The value of the objective function at x_k.
%
% iter       The number of iterations used to solve the problem.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Jan 8, 2003.  Last modified Feb 6, 2003.

[Inform, x_k, Obj, Iter] = ...
    miqpBBd(full(A), bl, bu, IntVars, Priority, ...
	    'HxFunc', mlp, kmax, nStackMax, optPar, ...
	    PriLev, PrintFile, Prob, moremem);

% MODIFICATION LOG
%
% 030109 ango Wrote file
% 030114 ango Revised comments (iprint/PriLev)
% 030206 hkh  Expand optPar to 20

% filterSQP NLP Matlab Solver
% -----------------------------------------------------------------------
%
% filterSQP solves the following constrained nonlinear programming
% problem (NLP):
%
%
% minimize f(x) subject to
%
%         (   x  )
%  bl <=  (  Ax  )  <= bu
%         ( c(x) )
%
%  where
%
%  f(x) is a nonlinear function of n variables.
%  c(x) is a nnCon vector of nonlinear constraint functions.
%  A    is a nnLin * n dense or sparse matrix with linear
%       constraint coefficient.
%
%  bl,bu have dimension n+m (m=nnCon+nnLin) and the elements are
%        ordered [bounds ; nonlinear ; linear]
%
%  The full input matrix A has three parts:
%
%  A = [g ConsPattern' A']
%
%  where g is a vector of length n, values irrelevant, ConsPattern is
%  the 0-1 pattern of the nonlinear constraint gradients and A is the
%  linear constraint coefficient matrix.
%
% -------------------------------------------------------------------
%   IF RUNNING TOMLAB:
%
%  Use GUI, driver routine tomRun or call filterSQPTL.m
%
% -------------------------------------------------------------------
% function [ifail,x_k,f_k,c_k,v_k,lws,istat,rstat] = filSQPs(...
%     A, bl, bu, nnCon, x_0, Scale, scmode, fLow, MaxIter,...
%     rho, mlp, kmax, maxf, WarmStart, lws, istat,...
%     PriLev, pname, optPar, Prob, moremem);
%
%
% The sparse version MEX is filSQPs, the dense is filSQPd
%
% -------------------------------------------------------------------
% INPUT:
%
% A          Gradient matrix [g ConsPattern' A'] (sparse or dense).
%
% bl         Lower bounds on (x,c(x),Ax).
% bu         Upper bounds on (x,c(x),Ax).
%
% nnCon      Number of nonlinear constraints (i.e. length(c(x)).
%
% x_0        Initial x vector (if empty set as 0) .
%
% Scale      n+m vector scale factors for variables and constraints
%            (same ordering as bl,bu).
%
% scmode     Scale mode:
%
%        0 - unit variable and constraint scaling (Scale can be
%            set empty)
%
%        1 - User provided scale factors for variables. Scale
%            must be of length n
%
%        2 - Unit variable scaling, user provided constraint
%            scaling. Scale must be of length n+m, but only the
%            last m elements are used.
%
%        3-  User provided variable AND constraint scaling. Scale
%            must be of length n+m (n+nnCon+nnLin)
%
%
% fLow       A lower bound on the objective function value.
%
% MaxIter    Maximum number of iterations.
%
% rho        Initial trust-region radius.
%
% mlp        Maximum level parameter for resolving degeneracy in BQPD.
%
% kmax       Maximum size of null-space (at most n).
%
% maxf       Maximum size of the filter.
%
% WarmStart  Set to 1 to restart the solver. If a warmstart is
%            requested, the input parameters lws and istat must be
%            provided. Also, n and m (the number of variables and
%            constraints) may not change.
%
% lws        Used only when doing a warmstart. This must be the
%            lws vector returned by the previous call to filterSQP.
%            Otherwise, set to empty.
%
% lam        Multipliers, n+m values required for warmstarts.
%            If wrong length, zeros are set in the solver.
%
% istat      Used only when doing a warmstart. Must be the first
%            element of the istat vector returned by the previous
%            call to filterSQP. Otherwise, set to empty.
%
% PriLev     Print level. Also see input parameter pname.
%        0 - Silent, except for minor output into <pname>.out.
%        1 - One line per iteration
%        2 - Scalar information printed
%        3 - Scalar and vector information printed
%
% pname      Problem name, at most 10 characters. The output files
%            are named <pname>.sum and <pname>.out
%
%    optPar     Vector of max length 20 with optimization parameters:
%               If any element is -999, default value is assigned.
%
%            Elements 2-8 are BQPD parameters, 1,9-11,19-20 for filterSQP.
%
%    optPar(1):  iprint  0     Print level in filterSQP
%    optPar(2):  tol     1E-10 Relative tolerance for BQPD subsolver
%    optPar(3):  emin    1.0   1=Use cscale in BQPD, 0=no scaling
%    optPar(4):  sgnf    5E-4  Max rel error in two numbers equal in
%                              exact arithmetic (BQPD)
%    optPar(5):  nrep    2     Max number of refinement steps (BQPD)
%    optPar(6):  npiv    3     No repeat if no more than npiv steps were taken
%    optPar(7):  nres    2     Max number of restarts if unsuccessful
%    optPar(8):  nfreq   500   The max interval between refactorizations
%
%    optPar(9):  NLP_eps 1E-6  NLP subproblem tolerance
%    optPar(10): ubd     1E-2  Upper bound on constraint violation used
%                              in the filter
%    optPar(11): tt      0.125 Parameter related to ubd. The actual
%                              upper bound is defined by the maximum of
%                              ubd and tt multiplied by the initial
%                              constraint violation
%
%    optPar(19): infty   1E20  A large value representing infinity
%
%    optPar(20): Nonlin  0     If 1, skip linear feasibility tests
%                              filterSQP treating all constraints as nonlinear
%
%
% Prob       The Tomlab problem definition structure.
%
% moremem   Scalar or 2x1-vector with workspace increase.
%           If <0, use default strategy.
%           If scalar, use same increase for both real and integer workspaces.
%           If vector, first element is for real workspace, second for integer.
%
%
% ---------------------------------------------------------------------
% OUTPUTS:
%
% Inform     Exit flag indicating success or failure:
%
%        0 - Solution found
%        1 - Unbounded: feasible point x with f(x)<=fmin found
%        2 - Linear constraints are infeasible
%        3 - (Locally) nonlinear infeasible, optimal solution to
%            feasibility problem found
%        4 - Terminated at point with h(x)<=eps but QP infeasible
%        5 - Terminated with rho<=eps
%        6 - Terminated due to too many iterations
%        7 - Crash in user routine could not be resolved
%        8 - Unexpected failure in QP solver
%        9 - Not enough real workspace
%       10 - Not enough integer workspace
%
% x_k        Solution vector.
% f_k        Function value at optimum x_k.
% c_k        Nonlinear constraints vector at optimum.
% v_k        Lagrange multipliers vector (bounds, nonlinear,
%            linear).
%
% lws        Integer vector (used as input when doing warmstarts).
%
% istat      Solution statistics, integer values. First element is
%            required as input if doing a warmstart.
%
%   istat(1)  Dimension of nullspace at solution
%   istat(2)  Number of iterations
%   istat(3)  Number of feasibility iterations
%   istat(4)  Number of objective evaluations
%   istat(5)  Number of constraint evaluations
%   istat(6)  Number of gradient evaluations
%   istat(7)  Number of Hessian evaluations
%   istat(8)  Number of QPs with mode<=2
%   istat(9)  Number of QPs with mode>=4
%   istat(10) Total number of QP pivots
%   istat(11) Number of SOC steps
%   istat(12) Maximum size of filter
%   istat(13) Maximum size of Phase 1 filter
%   istat(14) Number of QP crashes
%
%
% rstat      Solution statistics, real values.
%
%    rstat(1) l_2 norm of KT residual
%    rstat(2)
%    rstat(3) Largest modulus multiplier
%    rstat(4) l_inf norm of final step
%    rstat(5) Final constraint violation h(x)
%    rstat(6)
%    rstat(7)

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written Dec 13, 2002.  Last modified Dec 10, 2006.

function [Inform,x_k,f_k,c_k,v_k,lws,istat,rstat] = filSQPs(...
    A, bl, bu, nnCon, x_0, Scale, scmode, fLow, MaxIter,...
    rho, mlp, kmax, maxf, WarmStart, lws, lam, istat,...
    PriLev, pname, optPar, Prob, moremem)

help filterSQP;
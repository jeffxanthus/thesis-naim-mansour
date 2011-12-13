% BQPD QP Solver
% -----------------------------------------------------------------------
%
%   bqpd solves the following quadratic programming problem (QP):
%
%       minimize  0.5 x' H x + c'x   subject to:
%          x                                      (  x )
%                                          bl <=  ( Ax )  <= bu
%   where
%
%   H is an n x n sparse or dense matrix. Empty if LP problem.
%   A is an mA x n sparse matrix (linear constraints)
%   c is an n x 1 dense vector of linear objective coefficients.
%
%   The full input matrix A has two parts A = [c, A'];
%
%   bl, bu have dimension m=n+mA
%
%   TOMLAB: bl, bu and other input are generated in bqpdTL.m from Prob input.
%
% -----------------------------------------------------------------------
%   IF RUNNING TOMLAB:
%
%   Use driver routine tomRun or call bqpdTL.m.
% -----------------------------------------------------------------------
%
% function ...
%    [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpds(A, x_0, ...
%       bl, bu, H, fLowBnd,mlp,mode,kmax,PriLev,PrintFile, ...
%       k,ls,e,peq,lp,optPar,Prob,moremem);
%
% The sparse version MEX is bqpds, the dense is bqpdd
%
% ------------------------------------------------------------------------
% INPUT:
%
% A         Constraint matrix, n x m+1 (SPARSE).
% x_0       Initial x vector (if empty set as 0).
% bl        Lower bounds on (x,Ax).
% bu        Upper bounds on (x,Ax).
% H         Quadratic matrix, n x n, SPARSE  or DENSE, empty if LP problem.
%           If H is a string, H should be the name of a function routine, e.g
%           if H = 'HxComp' then the function routine
%
%           function Hx = HxComp(x, nState, Prob)
%
%           should compute H * x. The user must define this routine
%           nState == 1 if calling for the first time, otherwise 0.
%           Third argument, the Prob structure, should only be used if
%           calling BQPD with the additional input parameter Prob, see below.
%
%           Tomlab implements this callback to the predefined Matlab function
%           HxFunc.m, using the call if Prob.DUNDEE.callback == 1.
%
% fLowBnd   Lower bound on optimal f(x).
% mlp       Maximum number of levels of recursion.
% mode      Mode of operation, default set as 2*Prob.WarmStart.
% kmax      Max dimension of reduced space (k), default n, set as 0 if LP.
% PriLev    Print Level.
%           (0 = off, 1 = summary, 2 = scalar information, 3 = verbose).
% PrintFile Name of the Print file. Fortran unit 9 is used.
%           Name includes the path, maximal number of characters = 500.
%           Output is written on file bqpd.txt, if not given.
% To make bqpd to not open and not write anything to file: Set PriLev = 0.
%
% For Warm Start:
% k         Dimension of the reduced space (Warm Start).
% ls        Indices of active constraints, first n-k used for warm start.
% e         Steepest-edge normalization coefficients (Warm Start).
% peq       Pointer to the end of equality constraint indices in ls (Warm Start).
% lp        List of pointers to recursion information in ls (Warm Start).
%
% optPar    Vector of optimization parameters. If -999, set to default.
%           Length from 0 to 20 allowed.
%
% optPar(1): iprint 0    Print level in BQPD
% optPar(2): tol   1E-10 Relative accuracy in solution
% optPar(3): emin  1.0   1.0 Use cscale (constraint scaling) 0.0 no scaling
% optPar(4): sgnf  5E-4  Max rel error in two numbers equal in exact arithmetic
% optPar(5): nrep  2     Max number of refinement steps
% optPar(6): npiv  3     No repeat if no more than npiv steps were taken
% optPar(7): nres  2     Max number of restarts if unsuccessful
% optPar(8): nfreq 500   The max interval between refactorizations
%
% Prob      Sending the Prob structure is optional, only of use if sending
%           H as a function string, see input H.
%
% moremem   Scalar or 2x1-vector with workspace increase.
%           If <0, use default strategy.
%           If scalar, use same increase for both real and integer workspaces.
%           If vector, first element is for real workspace, second for integer.
%
% ----------------------------------------------------------------------
%
% OUTPUT:
%    [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpds(A, x_0, ...
%
% Inform    Result of BQPD run.
%           0 = Optimal solution found
% x_k       Solution vector with n decision variable values.
% Obj       Objective function value at optimum. If infeasible, the sum of
%           infeasibilities.
% g         Gradient at solution.
% Iter      Number of iterations.
%
% For Warm Start:
%
% k         Dimension of the reduced space (Warm Start).
% ls        Indices of active constraints, first n-k used for warm start.
% e         Steepest-edge normalization coefficients (Warm Start).
% peq       Pointer to the end of equality constraint indices in ls (Warm Start).
% lp        List of pointers to recursion information in ls (Warm Start).
%
% v_k       Lagrange parameters.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 26, 2002.   Last modified Jan 3, 2004.

function ... 
   [Inform, x_k, Obj, g, Iter,  k, ls, e, peq, lp, v_k] = bqpds(A, x_0, ...
      bl, bu, H, fLowBnd,mlp,mode,kmax,PriLev,PrintFile, ...
      k,ls,e,peq,lp,optPar,Prob);

help bqpd;
% minlpBB MINLP Matlab Solver
% -----------------------------------------------------------------------
%
% minlpBB solves the following constrained mixed-integer nonlinear
% programming problem (MINLP):
%
%
% minimize f(x) subject to
%
%         (   x  )
%  bl <=  (  Ax  )  <= bu
%         ( c(x) )
%
%  x_i integer valued for i in IntVars
%
%  where
%
%  f(x) is a nonlinear function of n variables.
%  c(x) is a nnCon vector of nonlinear constraint functions.
%  A    is a nnLin * n dense or sparse matrix with linear
%       constraint coefficients.
%
%  bl,bu have dimension n+m (m=nnCon+nnLin) and the elements are
%        ordered [bounds ; nonlinear ; linear]
%
%  The full input matrix A has three parts:
%
%  A = [g ConsPattern' A']
%
%  where g is a vector of length n, values are irrelevant, ConsPattern
%  is the 0-1 pattern of the nonlinear constraint gradients and A is
%  the linear constraint coefficient matrix.
%
% -------------------------------------------------------------------
%  IF RUNNING TOMLAB:
%
%  Use driver routine tomRun or call minlpBBTL.m
%
% -------------------------------------------------------------------
%
% function ...
%   [Inform,x_k,f_k,c_k,v_k,iter_nlp] = minlpBBs(...
%      A,bl,bu,nnCon,x_0,Scale,scmode,fLow,MaxIter,rho,...
%      kmax,mlp,maxf,PriLev,pname,optPar,Prob,IntVars,Priority,smax,...
%      setbeg,setcolidx,setrefrow,setprio,moremem);
%
% The sparse version MEX is minlpBBs, the dense is minlpBBd
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
% x_0        Initial x vector (if empty set as 0).
%
% Scale      n+m vector scale factors for variables and constraints
%            (ordered as bl,bu).
% scmode     Scale mode:
%
%        0 - unit variable and constraint scaling (Scale can be
%            set empty).
%
%        1 - User provided scale factors for variables. Scale
%            must be of at least length n.
%
%        2 - Unit variable scaling, user provided constraint
%            scaling. Scale must be of length n+m but only the last
%            m elements are used.
%
%        3-  User provided variable AND constraint scaling. Scale
%            must be of length n+m (n+nnCon+nnLin).
%
% MaxIter    Maximum number of NLP iterations for each node in the
%            branch-and-bound search tree.
%
% rho        Initial trust-region radius.
%
% kmax       Maximum size of null-space (at most n).
%
% maxf       Maximum size of the filter.
%
% mlp        Maximum level parameter for resolving degeneracy in BQPD
%            QP subsolver.
%
% PriLev     Print level. Also see input parameter pname.
%        0 - Silent, except for minor output into <pname>.out.
%        1 - Summary information written to <pname.out>
%        2 - Extensive information written to <pname.out>
%       >3 - Same as 2, but NLP subsolver filterSQP is called with
%            print level PriLev-3. Write to the same files.
%
% pname      Problem name, at most 10 characters. The out files
%            are named <pname>.sum and <pname>.out.
%
% optPar     Vector of max length 20 with optimization parameters.
%            If any element is -999, default value is assigned.
%            The elements used by minlpBB are:
%
%            Elements 1-7 are BQPD parameters, 8-10 for filterSQP,
%            11-13,17,20 for minlpBB
%
% optPar(1):  iprint   0     Print level in minlpBB
%
% optPar(2):  tol      1E-10 Relative tolerance for BQPD subsolver
% optPar(3):  emin     1.0   1=Use cscale in BQPD, 0=no scaling
% optPar(4):  sgnf     5E-4  Max rel error in two numbers equal in
%                            exact arithmetic (BQPD)
% optPar(5):  nrep     2     Max number of refinement steps (BQPD)
% optPar(6):  npiv     3     No repeat if no more than npiv steps were taken
% optPar(7):  nres     2     Max number of restarts if unsuccessful
% optPar(8):  nfreq    500   The max interval between refactorizations
%
% optPar(9):  NLP_eps  1E-6  NLP subproblem tolerance
% optPar(10): ubd      1E-2  Upper bound on constraint violation used
%                            in the filter
% optPar(11): tt       0.125 Parameter related to ubd. The actual
%                            upper bound is defined by the maximum of
%                            ubd and tt multiplied by the initial
%                            constraint violation
%
% optPar(12): epsilon  1D-6  Tolerance for x-value tests
% optPar(13): MIopttol 1E-4  Tolerance for function value tests
%
% optPar(17): branchtype     Branch strategy, currently only 1 strategy
%
% optPar(19): infty    1E20  A large value representing infinity
%
% optPar(20): Nonlin   0     If 1, skip linear feasibility tests
%                            minlpBB treating all constraints as nonlinear
%
%
% Prob       The Tomlab problem definition structure. The names of
%            the m-files calculating function and constraint values,
%            gradients and second derivatives must be present in Prob.
%
% IntVars    Vector of indices of the integer variables.
%
% Priority   Vector of priorities for the integer variables. Higher
%            values imply higher priority (NOTE: priorities must be
%            integer values, any fractional parts are truncated).
%
% smax       Maximum size of the stack storing information during
%            the tree search.
%
%  --- SOS1 set parameters, setbeg, setcolidx, setrefrow, setprio ---
%
% setbeg     Vector of pointers to start of each SOS1 set in
%            setcolidx and setrefrow. One extra element must be
%            added at the end, set to length(setbeg)+1.
%
% setcolidx  The indices of the variables in each SOS1 set.
%            Different sets are separated by the pointers in setbeg.
%
% setrefrow  Reference rows with ordering information for the
%            variables in.
%
% setprio    Priorities among SOS1 sets (Integer values used).
%
%
% moremem    Scalar or 2x1-vector with workspace increase.
%            If <0, use default strategy.
%            If scalar, use same increase for both real and integer workspaces.
%            If vector, first element is for real workspace, second for integer.
%
% fLow       Lower bound for QP subproblems.
%
% lam        Vector of multipliers. Loaded if length is n+m,
%            otherwise ignored.
%
% ---------------------------------------------------------------------
% OUTPUT:
%
% ifail      minlpBB information parameter.
%
% x_k        Solution vector.
% f_k        Function value at optimum.
% c_k        Constraint residuals.
% v_k        Lagrangian multipliers (for bounds + dual solution vector).
%
% iter_nlp   Number of iterations.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Dec 13, 2002.  Last modified Dec 10, 2006.

function [ifail,x_k,f_k,c_k,v_k,iter_nlp] = minlpBB(...
    A,bl,bu,nnCon,x_0,Scale,scmode,MaxIter,rho,...
    kmax,mlp,maxf,PriLev,pname,optPar,Prob,IntVars,Priority,...
    smax,setbeg,setcolidx,setrefrow,setprio,fLow,lam)
  
help minlpBB;

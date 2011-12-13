% NLSSOL NLLS Matlab Solver
% -----------------------------------------------------------------------
%
% nlssol  solves the constrained nonlinear least-squares problem
%
%         minimize           (1/2)(y - r(x))'(y - r(x))
%
%                                 (    x  )
%         subject to    bl  .le.  (  A*x  )  .le.  bu
%                                 (  c(x) )
%
% where y is a given constant m-vector, r(x) is an m-vector of
% smooth functions, A is a constant matrix and c(x) is a vector
% of smooth nonlinear functions. The feasible region is defined
% by a mixture of linear and nonlinear equality or inequality
% constraints on x.
%
% The dimensions of the problem are:
%
% m        the number of data points (the dimension of y),
%
% n        the number of variables (dimension of  x),
%
% nclin    the number of linear constraints (rows of the matrix  A),
%
% ncnln    the number of nonlinear constraints (dimension of  c(x)).
%
% bl, bu have dimension nb = n+nclin+ncnln
%
% nlssol uses a sequential quadratic programming algorithm, with a
% positive-definite quasi-Newton approximation to the transformed
% Hessian  Q'HQ  of the Lagrangian function (which will be stored in
% the array R).
% -----------------------------------------------------------------------
%
% IF RUNNING TOMLAB:
%
% bl, bu and other input are generated in nlssolTL.m from Prob input.
% Define nonlinear residual and Jacobian and constraint and gradients for the
% constraints according to TOMLAB standard.
% Use driver routine tomRun, or nlssolTL, but do not call nlssol directly.
%
% IF NOT RUNNING TOMLAB:
%
% The file 'funfdf.m' must be defined and contain:
% function [r,J]=funfdf(x,Prob,mode,nstate) to compute the residual vector r
% and the Jacobian matrix J at the point x.
%
% The file 'funcdc.m' must be defined and contain:
% function [c,dc]=funcdc(x,Prob,mode,nstate) to compute the nonlinear
% constraint value c and the constraint Jacobian dc for
% the nonlinear constraints at the point x.
%
% See comments below for the INPUT variable Prob.
%
% -----------------------------------------------------------------------
%
% function ...
% [x, Inform, iState, cLamda, iwCount, fObj, gObj, r, J, fCon, gCon, H] = ...
%     nlssol( A, bl, bu, x, Prob, y, optPar, ...
%            Warm, H, iState,  cLamda, ...
%            SpecsFile, PrintFile, SummFile, ...
%            PriLev, ProbName );
% ------------------------------------------------------------------------
% INPUT:
%
% A         Constraint matrix, nb x n (DENSE).
% bl        Lower bounds on (x,Ax,c(x)).
% bu        Upper bounds on (x,Ax,c(x)).
% x         Initial x vector (n x 1).
%           If Warm start x must correspond to values in H and iState.
% Prob      Must be a structure. No check is made in the MEX interface.
%
%           If TOMLAB calls nlssol, then Prob is the standard
%           TOMLAB problem structure, otherwise the user should set:
%
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
%           Two user written routines must be written:
%
%           funfdf, name stored in Prob.FUNCS.rJ, with syntax
%                [mode, r, J]     = funfdf(x, Prob, mode, nstate)
%
%           funcdc, name stored in Prob.FUNCS.cdc, with syntax
%                [mode, c, dc]   = funcdc(x, Prob, mode, nstate)
%
%           NLSOL is calling the TOMLAB routines ls_rJ.m and nlp_cdc.m
%           in the callback, and they call funfdf and funcdc, respectively.
%
%           If these fields in Prob are empty (Prob.FUNCS.rJ, Prob.FUNCS.cdc),
%           the TOMLAB callback routines calls the usual function routines.
%           Then the Prob struct should be normally defined, and
%           the fields Prob.FUNCS.r, Prob.FUNCS.J, Prob.FUNCS.c, Prob.FUNCS.dc
%           be set in the normal way (e.g. by the routine tomFiles.m, or one
%           of the Assign-routines like clsAssign.m).
%
%           If the mode parameter is 0, funfdf should return r, otherwise
%           both r and the Jacobian matrix J.
%           If the mode parameter is 0, funcdc should return c, otherwise
%           both c and dc. Note that each row in dc corresponds to a
%           constraint, and that dc should be a dense matrix.
%           If the matrix dc is returned as a sparse Matlab matrix, nlp_cdc
%           will do full(dc) to get a dense matrix.
%
%           The user could also write his own versions of the routines
%           ls_rJ.m and nlp_cdc.m and put them before in the path.
% y         Data vector of length. Number of residuals == length(y) == m.
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           Set empty if only using default parameters.
% Warm      Flag for Warm start (==1) or Cold Start (==0 ), i.e. normal mode.
%           If 'Warm Start', iState, cLamda, H and x must supply correct values
% H         Hessian matrix, only accessed if Warm start.
%           Must select: Hessian = Yes in order to do a warm start.
%           Hessian Yes is equivalent to optPar(50) = 1 (default is 0).
% iState    Working set (if Warm start) (nb = n+nclin+ncnln) x 1.
%           If length(iState) < nb, setting iState(1:nb)=0;
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
% cLamda    Lagrangian multipliers for the n + nclin + ncnln  constraints.
%           If Warm start, cLamda(n+nclin+1:n+nclin+ncnln), the nonlinear
%           Lagrange multipliers, must correspond to values in iState.
% SpecsFile Name of the OPTIONS File, see TOMLAB /SOL guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number
%           of characters = 500.
% PriLev    Print level in the nlssol MEX-interface.
%           = 0  Silent
%           = 10 Dimensions are printed
%           if isempty(PriLev), set as 0.
% ProbName  Name of the problem. <=100 characters are used in the MEX interface
%           (Not used by MEX-interface).
%
% Implicit input:
% m         Number of residuals, length(y).
% nb        length(bl).
% n         Number of nonlinear objective variables, length(x).
% nclin     Number of linear constraints size(A,1).
% ncnln     Number of nonlinear constraints, nb - length(x)-size(A,1).
%
% DESCRIPTION of optPar vector:
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by NLSSOL,
% if not the value is altered in the SPECS file (input SpecsFile). Refer to 
% nlssolTL.m for information about optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT:
% x         Solution vector (n by 1) with n decision variable values.
% Inform    Result of NLSSOL run (See NPSOL User's Guide).
%           0 = Optimal solution found
%           1 = Point satisfies Kuhn-Tucker conditions, but sequence
%               of iterates has not yet converged
%           2 = No feasible point found for given Linear feasibility tolerance
%           3 = No feasible point found for given nonlinear constraints
%           4 = Maximal number of iterations reached
%           6 = No sufficient decrease in merit function during line search
%           7 = Large errors in the derivatives
%           9 = An input parameter is invalid
% iState    Status of working set, se input description of iState.
% cLamda    Lagrangian multipliers (dual solution vector) (nb x 1 vector).
% iwCount   Number of iterations (iwCount(1)), function evaluations
%           (iwCount(2)) and constraint evaluations (iwCount(3)).
% fObj      Objective function value at optimum.
% gObj      Gradient of the objective, n x 1.
% r         Residual vector m x 1.
% J         Jacobian matrix m x n.
% fCon      Nonlinear constraint vector, ncnln x 1.
% gCon      Gradient matrix of the nonlinear constraint vector, ncnln x n.
% H         Cholesky factor of Hessian approximation.
%           Hessian no  - reordered variables.
%           Hessian yes - natural order of variables, used for Warm start.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2011 by Tomlab Optimization Inc., Sweden. $Release: 7.8.0$
% Written Sept 17, 2000.     Last modified July 24, 2011.

function ...
    [x, Inform, iState, cLamda, iwCount, fObj, gObj, r, J, fCon, gCon, H] = ...
    nlssol( A, bl, bu, x, Prob, y, optPar, ...
    Warm, H, iState, cLamda, ...
    SpecsFile, PrintFile, SummFile, ...
    PriLev, ProbName );

help nlssol;

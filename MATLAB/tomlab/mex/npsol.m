% NPSOL NLP Matlab Solver
% -----------------------------------------------------------------------
%
% npsol solves the following nonlinear programming problem (NLP):
%
%       minimize  f(x)            subject to:
%          x                               (  x    )
%                                   bl <=  (  A x  )  <= bu
%                                          (  c(x) )
% where
%
% A is an nclin x n dense matrix (linear constraints).
% f(x) has n nonlinear variables.
% c(x) is an ncnln x 1 vector of nonlinear constraints.
% The Jacobian of d/dx c(x) is a ncnln x n dense matrix.
%
% bl, bu have dimension m = n+nclin+ncnln.
%
% -----------------------------------------------------------------------
%
% IF RUNNING TOMLAB:
%
% bl, bu and other input are generated in npsolTL.m from Prob input.
% Define nonlinear functions and gradients for the
% function and constraints according to TOMLAB standard.
% Use driver routine tomRun or npsolTL, but do not call npsol directly
%
% IF NOT RUNNING TOMLAB:
%
% The file 'funfdf.m' must be defined and contain:
% function [mode,f,g]=funfdf(x,Prob,mode,nstate) to compute the objective
% function f and the gradient g at the point x.
%
% The file 'funcdc.m' must be defined and contain:
% function [mode,c,dc]=funcdc(x,Prob,mode,nstate) to compute the nonlinear
% constraint value c and the constraint Jacobian dc for
% the nonlinear constraints at the point x.
%
% See comments below for the INPUT variable Prob
%
% -------------------------------------------------------------------------
%
% function [x, Inform, iState, cLamda, iwCount, fObj, gObj, fCon, gCon, H] = ...
%     npsol( A, bl, bu, x, Prob, optPar, ...
%            Warm, H, iState, cLamda, ...
%            SpecsFile, PrintFile, SummFile, ...
%            PriLev, ProbName );
%
% ------------------------------------------------------------------------
% INPUT:
%
% A         Constraint matrix, m x n (DENSE).
% bl        Lower bounds on (x,Ax,c(x)).
% bu        Upper bounds on (x,Ax,c(x)).
% x         Initial x vector (n x 1).
%           If Warm start x must correspond to values in H and iState.
% Prob      Must be a structure. No check is made in the MEX interface.
%
%           If TOMLAB calls npsol, then Prob is the standard
%           TOMLAB problem structure, otherwise the user should set:
%
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
%           Two user written routines must be written:
%
%           funfdf, actual name stored in Prob.FUNCS.fg, with syntax
%                [mode, f, g]     = funfdf(x, Prob, mode, nstate)
%
%           funcdc, actual name stored in Prob.FUNCS.cdc, with syntax
%                [mode, c, dc]   = funcdc(x, Prob, mode, nstate)
%
%           NPSOL is calling the TOMLAB routines nlp_fg.m and nlp_cdc.m
%           in the callback, and they call funfdf and funcdc, respectively.
%
%           If these fields in Prob are empty (Prob.FUNCS.fg, Prob.FUNCS.cdc),
%           the TOMLAB callback routines calls the usual function routines.
%           Then the Prob struct should be normally defined, and
%           the fields Prob.FUNCS.f, Prob.FUNCS.g, Prob.FUNCS.c, Prob.FUNCS.dc
%           be set in the normal way (e.g. by the routine tomFiles.m, or one
%           of the Assign-routines like conAssign.m).
%
%           If the mode parameter is 0, funfdf should return f, otherwise
%           both f and the gradient vector g.
%           If the mode parameter is 0, funcdc should return c, otherwise
%           both c and dc. Note that each row in dc corresponds to a
%           constraint, and that dc should be a dense matrix.
%           If the matrix dc is returned as a sparse Matlab matrix, nlp_cdc
%           will do full(dc) to get a dense matrix.
%
%           The user could also write his own versions of the routines
%           nlp_fg.m and nlp_cdc.m and put them before in the path.
%
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           Set empty if only using default parameters
% Warm      Flag for Warm start (==1) or Cold Start (==0 ), i.e. normal mode.
%           If 'Warm Start', iState, cLamda, H and x must supply correct values.
% H         Hessian matrix, only accessed if Warm start.
%           Must select: Hessian = Yes in order to do a warm start.
%           Hessian Yes is equivalent to optPar(50) = 1 (default is 0).
% iState    Working set (if Warm start) (nb = n+nclin+ncnln) x 1 (DENSE).
%           If length(iState) < nb, setting iState(1:nb)=0;
% iState(i)=0: Corresponding constraint not in the initial working set.
% iState(i)=1: Inequality constraint at its lower bound in working set.
% iState(i)=2: Inequality constraint at its upper bound in working set.
% iState(i)=3: Equality constraint in the initial working set, bl(i)==bu(i).
% cLamda    Lagrangian multipliers for the n + nclin + ncnln constraints.
%           If Warm start, cLamda(n+nclin+1:n+nclin+ncnln), the nonlinear
%           Lagrange multipliers, must correspond to values in iState.
% SpecsFile Name of the input parameter file, see TOMLAB /SOL guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number
%           of characters = 500.
% PriLev    Print level in the npsol MEX-interface.
%           = 0  Silent
%           = 10 Dimensions are printed
% ProbName  Name of the problem. <=100 characters are used in the MEX interface
%           (Not used by MEX-interface).
%
% Implicit input:
% nb        length(bl) (== m).
% n         Number of nonlinear objective variables, length(x).
% nclin     Number of linear constraints, size(A,1).
% ncnln     Number of nonlinear constraints, nb - length(x)-size(A,1).
%
% DESCRIPTION of optPar vector:
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by NPSOL,
% if not the value is altered in the SPECS file (input SpecsFile). Refer to 
% npsolTL.m for information about optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT:
% x         Solution vector (n by 1) with n decision variable values.
% Inform    Result of NPSOL run (See NPSOL User's Guide).
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
% cLamda    Lagrangian multipliers (dual solution vector) (m x 1 vector).
% iwCount   Number of iterations (iwCount(1)),
%           function (iwCount(2)) and constraint (iwCount(3)) calls.
% fObj      Objective function value at optimum.
% gObj      Gradient of the nonlinear objective.
% fCon      Nonlinear constraint vector.
% gCon      Gradient matrix of the nonlinear constraint vector.
% H         Cholesky factor of Hessian approximation.
%           Hessian no  - reordered variables.
%           Hessian yes - natural order of variables, used for Warm start.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization AB, Sweden. $Release: 5.5.0$
% Written Sept 15, 2000.     Last modified Aug 14, 2006.

function [x, Inform, iState, cLamda, iwCount, fObj, gObj, fCon, gCon, H] = ...
    npsol( A, bl, bu, x, Prob, optPar, ...
    Warm, H, iState, cLamda, ...
    SpecsFile, PrintFile, SummFile, ...
    PriLev, ProbName );

help npsol;

fprintf('\n\n\n')
fprintf('------------------------------------------------------------\n')
fprintf('NOTE!!! If this help executes, it means that no npsol.dll \n')
fprintf('is in the PATH to tomlab. Your installation is wrong!!!\n')
fprintf('------------------------------------------------------------\n')
fprintf('\n\n\n')
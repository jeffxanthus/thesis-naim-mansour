% MINOS NLP/QP/LP Matlab Solver
% -----------------------------------------------------------------------
%
%   minos either solves the following nonlinear programming problem (NLP):
%
%        minimize  f(x) + c'x    subject to:
%           x                           (  x    )
%                                bl <=  (  g(x) )  <= bu
%                                       (  Ax   )
%   where
%
%   A is an mA x n sparse matrix (linear constraints)
%   c is a n x 1 vector of linear objective coefficients.
%   f(x) has nnObj nonlinear variables
%   g(x) has dimension nnCon.
%   The Jacobian is d/dx g(x) has nnJac variables, and ne nonzero elements
%
%   or the quadratic programming problem (QP):
%
%        minimize  0.5 * x' * H * x  + c'x    subject to:
%           x                           (  x    )
%                                bl <=  (  g(x) )  <= bu
%                                       (  Ax   )
%
%   Full input matrix A has three parts A = [d/dx g(x); A; c'];
%   The position of the row c' is iObj. iObj=0 means no linear part
%
%   bl, bu have dimension m=n+nnCon+mA+(iObj~=0)
%
%   If nnObj = nnJac = nnCon = 0, and g(x) is empty, a standard sparse (LP)
%   is solved in an efficient way.
%
%   Note that for a standard QP g(x) is empty. However, as the interface
%   is implemented, nonempty g(x) is handled the same way as for general NLP
%
% -----------------------------------------------------------------------
%
%   IF RUNNING TOMLAB:
%
%   bl and bu and other inputs are generated in minosTL for NLP problems (or
%   minosqpTL for QP problems, or minoslpTL for LP problems) from Prob input.
%   Define nonlinear functions, constraints and gradients for the
%   function and constraints according to TOMLAB standard.
%   Use driver routine tomRun, or call minosTL, minoslpTL or minosqpTL. 
%   Do not call minos.m directly.
%
%   IF NOT RUNNING TOMLAB:
%
%   The file 'funfdf.m' must be defined and contain:
%   function [mode,f,g]=funfdf(x,Prob,mode,nstate) to compute the objective
%   function f and the gradient g at the point x.
%
%   The file 'funcdc.m' must be defined and contain:
%   function [mode,c,dcS]=funcdc(x,Prob,mode,nstate) to compute the nonlinear
%   constraint value c and the constraint Jacobian dcS for
%   the nonlinear constraints at the point x.
%
%   NOTE: The matrix dcS MUST be a SPARSE MATLAB matrix.
%   Do dcS = sparse(dcS); after dcS has been computed.
%
%   See comments below for the INPUT variable Prob.
%
% ------------------------------------------------------------------------
%
% [hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
%      minos( H, A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
%             Warm, hs, xs, pi, nS, ...
%             SpecsFile, PrintFile, SummFile, ...
%             PriLev, ObjAdd, moremem, ProbName );
% ------------------------------------------------------------------------
% INPUT:  Must give at least 8 first parameters.
%
% H         Matrix n x n in a quadratic programming (QP) problem.
%           DENSE or SPARSE. Leave empty if LP, or NLP problem.
% A         Constraint matrix, m x n SPARSE (nonlinear, linear and objective)
%           m > 0 always!!! Define dummy constraint for unconstrained
%           problems.
% bl        Lower bounds on (x,g(x),Ax,c').
% bu        Upper bounds on (x,g(x),Ax,c').
%           NOTE! The bl and bu values for the last nonlinear constraint c
%           must have reverse signs and be put in each other places:
%           If     c_L   <=  c(x)   <= c_U, then bl = -c_U  and bu = -c_L.
%           This is because the bounds acts as the constraints on the
%           slack variables for the nonlinear constraints.
% nnCon     Number of nonlinear constraints.
% nnObj     Number of nonlinear objective variables.
% nnJac     Number of nonlinear Jacobian variables.
%
% Prob      Must be a structure. No check is made in the MEX interface.
%
%           If TOMLAB calls minos, then Prob is the standard
%           TOMLAB problem structure, otherwise the user should set:
%
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
%           If the problem is a LP or QP problem (H defined), the user
%           does not have to specify anything else in the structure.
%
%           For a general nonlinear objective or nonlinear constraints
%           names of two user written routines must be given:
%
%           funfdf, actual name stored in Prob.FUNCS.fg, with syntax
%                [mode, f, g]     = funfdf(x, Prob, mode, nstate)
%
%           funcdc, actual name stored in Prob.FUNCS.cdc, with syntax
%                [mode, c, dcS]   = funcdc(x, Prob, mode, nstate)
%
%           MINOS is calling the TOMLAB routines nlp_fg.m and nlp_cdcS.m
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
%           both c and dcS. Note that each row in dcS corresponds to a
%           constraint, and that dcS must be a SPARSE matrix.
%
%           The user could also write his own versions of the routines
%           nlp_fg.m and nlp_cdcS.m and put them before in the path.
%
% iObj      Says which row of A is a free row containing a linear objective
%           vector c. If there is no such vector, iObj = 0. Otherwise, this
%           row must come after any nonlinear rows, so that nnCon <= iObj <= m
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           If using only default options, set optPar as an empty matrix.
% Warm      Flag, if true: warm start. Default cold start (if empty)
%           If 'Warm Start' xs, nS and hs must be supplied with correct
%           values.
% hs        Basis status of variables + constraints (n+m x 1 vector)
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%           2=superbasic (between bounds), 3=basic (between bounds).
% xs        Initial vector, optionally including m slacks at the end.
%           If warm start, full xs must be supplied.
% pi        Lagrangian multipliers for the nnCon nonlinear constraints.
%           If empty, set as 0.
% nS        # of superbasics. Only used if calling again with a Warm Start.
% SpecsFile Name of the SPECS input parameter file, TOMLAB /MINOS guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number 
%           of characters = 500.
% PriLev    Printing level in the minos m-file and minos MEX-interface.
%           = 0  Silent
%           = 1  Summary information
%           = 2  More detailed information
% ObjAdd    Constant added to the objective for printing purposes, typically 0.
% moremem   Add extra memory for the sparse LU, might speed up the optimization.
%           1E6 is 10MB of memory. If empty, set as 0.
% ProbName  Name of the problem. <=100 characters are used in the MEX interface.
%           In the MINOS solver the first 8 characters are used in the printed
%           solution. Blank is OK.
%
% DESCRIPTION of optPar vector:
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by MINOS, unless the value is 
% altered in the SPECS file. Refer to minosTL.m for information about 
% optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT:
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%           2=superbasic (between bounds), 3=basic (between bounds).
% xs        Solution vector (n+m by 1) with n decision variable values
%           together with the m slack variables.
% pi        Lagrangian multipliers (dual solution vector) (m x 1 vector).
% rc        Reduced costs, a n+m vector. If nInf=0, last m == -pi
% Inform    Result of MINOS run.
%           0 = Optimal solution found
% nS        # of superbasics.
% nInf      Number of infeasibilities.
% sInf      Sum of infeasibilities.
% Obj       Objective function value at optimum.
% iwCount   Number of iterations (major and minor),
%           function and constraint calls.
% gObj      Gradient of the nonlinear objective.
% fCon      Nonlinear constraint vector.
% gCon      Gradient vector (non-zeros) of the nonlinear constraint vector.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written July 10, 2000.     Last modified Dec 24, 2006.
%# mex

function ...
[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     minos( H, A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, ProbName )

help minos;
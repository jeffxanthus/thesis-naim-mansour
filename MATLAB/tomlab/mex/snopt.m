% SNOPT NLP Matlab Solver
% -----------------------------------------------------------------------
%
%   snopt solves the following nonlinear programming problem (NLP):
%
%         minimize  f(x) + c'x      subject to:
%            x                             (  x    )
%                                   bl <=  (  g(x) )  <= bu
%                                          (  Ax   )  
%   where 
%
%   A is an mA x n sparse matrix (linear constraints).
%   c is an n x 1 vector of linear objective coefficients. 
%   f(x) has nnObj nonlinear variables.
%   g(x) is an nnCon x 1 vector of nonlinear constraints.
%   The Jacobian of d/dx g(x) has nnJac variables, and ne nonzero elements.
%
%   The full input matrix A has three parts A = [d/dx g(x); A; c'];
%   The position of the row c' is iObj. iObj=0 means no linear objective part.
%
%   bl, bu have dimension m=n+nnCon+mA+(iObj~=0).
%
% -----------------------------------------------------------------------
%
%   IF RUNNING TOMLAB:
%
%   bl, bu and other input are generated in snoptTL from Prob input.
%   Define nonlinear functions, constraints and gradients for the
%   function and constraints according to TOMLAB standard.
%   Use tomRun or snoptTL, but do not call snopt directly.
%
%   IF NOT RUNNING TOMLAB, and calling snopt directly:
%
%   The file 'funfdf.m' must be defined and contain:
%   function [mode, f,g]=funfdf(x,Prob,mode,nstate) to compute the objective
%   function f and the gradient g at the point x.
%
%   The file 'funcdc.m' must be defined and contain:
%   function [mode, c,dcS]=funcdc(x,Prob,mode,nstate) to compute the nonlinear
%   constraint value c and the constraint Jacobian dcS for
%   the nonlinear constraints at the point x.
%
%   Note that Matlab has dynamic sparse matrix handling and Fortran has
%   static handling. The returned vector of constraint gradient values
%   must always match the pattern of nonzeros as defined in the call to snopt.
%   One approach for a general solution to this is given in the Tomlab
%   callback routine nlp_cdcS.m, which calls the user defined 'funcdc' function
%
%   The fields Prob.P and Prob.ConsPattern must be set, see below.
%
%   See comments below for the INPUT variable Prob
%
% -------------------------------------------------------------------------
%
% function ...
% [hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
%      snopt( A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
%             Warm, hs, xs, pi, nS, ...
%             SpecsFile, PrintFile, SummFile, ...
%             PriLev, ObjAdd, moremem, ProbName );
% 
% ------------------------------------------------------------------------
% INPUT:  Must give at least 7 first parameters. 
%
% A         Constraint matrix, m x n SPARSE (A consists of nonlinear part, 
%           linear part and one row for the linear objective). m > 0 always.
% bl        Lower bounds on (x,g(x),Ax,c').
% bu        Upper bounds on (x,g(x),Ax,c').
% nnCon     Number of nonlinear constraints.
% nnObj     Number of nonlinear objective variables.
% nnJac     Number of nonlinear Jacobian variables.
%
% Prob      Must be a structure. No check is made in the MEX interface.
%
%           If TOMLAB calls snopt, then Prob is the standard
%           TOMLAB problem structure, otherwise the user should set:
%
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%
%           Prob.ConsPattern = []; 
%
%           or as the nonzero pattern for the constraint Jacobian as
%
%           Prob.ConsPattern = ConsPattern;
%
%           ConsPattern is a nnCon x n zero-one sparse or dense matrix, 
%           where 0 values indicate zeros in the constraint Jacobian and 
%           ones indicate values that might be non-zero. 
%
%           If the problem is a LP or QP problem (H defined), then the user 
%           does not have to specify anything more in the structure.
%
%           For a general nonlinear objective, or nonlinear constraints 
%           names of two user written routines must be given:
%
%           funfdf, actual name stored in Prob.FUNCS.fg, with syntax
%                [mode, f, g]     = funfdf(x, Prob, mode, nstate)
%
%           funcdc, actual name stored in Prob.FUNCS.cdc, with syntax 
%                [mode, c, dcS]   = funcdc(x, Prob, mode, nstate)
%
%           SNOPT is calling the TOMLAB routines nlp_fg.m and nlp_cdcS.m
%           in the callback, and they call funfdf and funcdc, respectively.
%
%           If these fields in Prob are empty (Prob.FUNCS.fg, Prob.FUNCS.cdc), 
%           the TOMLAB callback routines calls the usual function routines.
%           Then the Prob struct should be normally defined, and
%           the fields Prob.FUNCS.f, Prob.FUNCS.g, Prob.FUNCS.c, Prob.FUNCS.dc
%           be set in the normal way (e.g. by the routine tomFiles.m).
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
%           row must come after any nonlinear rows, so that nnCon <= iObj <= m.
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           If using only default options, set optPar as an empty matrix. 
% Warm      Flag, if true: warm start. Default cold start (if empty)
%           If 'Warm Start' xS, nS and hs must be supplied with correct values
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%                 2=superbasic (between bounds), 3=basic (between bounds)
% xs        Initial vector, optionally including m slacks at the end.
%           If warm start, full xs must be supplied.
% pi        Lagrangian multipliers for the nnCon nonlinear constraints.
%           If empty, set as 0.
% nS        # of superbasics. Only used if calling again with a Warm Start.
% SpecsFile Name of the SPECS input parameter file, see TOMLAB /SNOPT guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500.
% SummFile  Name of the Summary file. Name includes the path, maximal number 
%           of characters = 500.
% PriLev    Printing level in the snopt m-file and snopt MEX-interface.
%           = 0  Silent
%           = 1  Summary information
%           = 2  More detailed information
% ObjAdd    Constant added to the objective for printing purposes, typically 0.
% moremem   Add extra memory for the sparse LU, might speed up the optimization.
%           1E6 is 10MB of memory. If empty, set as 0.
% ProbName  Name of the problem. <=100 characters are used in the MEX interface
%           In the SNOPT solver the first 8 characters are used in the printed
%           solution. Blank is OK.
%
% DESCRIPTION of optPar vector:
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by SNOPT, if not the value is
% altered in the SPECS file (input SpecsFile). Refer to snoptTL.m for 
% information about optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT: 
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu)
%               2=superbasic (between bounds), 3=basic (between bounds).
% xs        Solution vector (n+m by 1) with n decision variable values 
%           together with the m slack variables.
% pi        Lagrangian multipliers (dual solution vector) (m x 1 vector).
% rc        Reduced costs, a n+m vector. If nInf=0, last m == pi.
% Inform    Result of SNOPT run.
%           0 = Optimal solution found
% nS        # of superbasics.
% nInf      Number of infeasibilities.
% sInf      Sum of infeasibilities.
% Obj       Objective function value at optimum.
% iwCount   Number of iterations minor (iwCount(1)) and major (iwCount(2)), 
%           function (iwCount(3:6)) and constraint (iwCount(7:10)) calls.
% gObj      Gradient of the nonlinear objective.
% fCon      Nonlinear constraint vector.
% gCon      Gradient vector (non-zeros) of the nonlinear constraint vector.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written July 4, 2000.   Last modified Dec 24, 2006.
%# mex

function ...
[hs, xs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount, gObj, fCon, gCon] = ...
     snopt( A, bl, bu, nnCon, nnObj, nnJac, Prob, iObj, optPar, ...
            Warm, hs, xs, pi, nS, ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, ProbName )

help snopt;
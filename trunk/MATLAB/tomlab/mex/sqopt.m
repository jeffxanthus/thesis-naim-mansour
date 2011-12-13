% SQOPT QP Matlab Solver
% -----------------------------------------------------------------------
%
%   sqopt solves the following quadratic programming problem (QP):
%
%       minimize  0.5 x' H x + c'x + d'x   subject to:
%          x                                      (  x )
%                                          bl <=  ( Ax )  <= bu
%   where 
%
%   H is an n x n sparse or dense matrix. Empty if LP problem.
%   A is an mA x n sparse matrix (linear constraints).
%   c is an n x 1 dense vector of linear objective coefficients. 
%   d is an n x 1 dense vector of linear objective coefficients. 
%
%   The full input matrix A has two parts A = [A; d'];
%   The position of the row d' is iObj. iObj=0 means no linear part in A.
%
%   bl, bu have dimension m=n+mA+(iObj~=0).
%
%   TOMLAB: bl, bu and other input are generated in sqoptTL.m from Prob input.
%
%   NOTE: There are two ways to give the linear objective: either explicit 
%   as vector c or as part of the sparse matrix A, as d (or both ways).
%
% -----------------------------------------------------------------------
%   IF RUNNING TOMLAB:
%
%   Use driver routine tomRun or sqoptTL, but do not call sqopt directly
% -----------------------------------------------------------------------
%
% function ... 
% [xs, hs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
%      sqopt( A, bl, bu, H, c, hElast, iObj, optPar, ...
%             Warm, hs, xs,  nS, ...
%             SpecsFile, PrintFile, SummFile, ...
%             ObjAdd, moremem, ProbName, Prob );
% 
% ------------------------------------------------------------------------
% INPUT:  
%
% A         Constraint matrix, m x n (SPARSE).
% bl        Lower bounds on (x,Ax,d').
% bu        Upper bounds on (x,Ax,d').
% H         Quadratic matrix, n x n, SPARSE  or DENSE, empty if LP problem
%           If H is a string, H should be the name of a function routine, e.g
%           if H = 'HxComp' then the function routine: 
%
%           function Hx = HxComp(x, nState, Prob)
%            
%           should compute H * x. The user must define this routine.
%           nState == 1 if calling for the first time, otherwise 0.
%           Third argument, the Prob structure, should only be used if
%           calling SQOPT with the additional input parameter Prob, see below
%
%           Tomlab implements this callback to the predefined Matlab function
%           HxFunc.m, using the call if Prob.SOL.callback == 1.
%
% c         Linear objective 
% hElast    defines which bounds are elastic in elastic mode. hElast(j):
%           0 = variable j cannot be infeasible.
%           1 = variable j can violate its lower bound.
%           2 = variable j can violate its upper bound.
%           3 = variable j can violate either its lower or upper bound.
% iObj      Says which row of A is a free row containing a linear objective
%           vector d. If there is no such vector, iObj = 0. 
% optPar    Vector with optimization parameters overriding defaults and the
%           optionally specified SPECS file.
%           Set empty if only using default parameters.
% Warm      Flag, if true: warm start. Default cold start (if empty)
%           If 'Warm Start' xS, nS and hs must be supplied with correct values
% hs        Basis status of variables + constraints (n+m x 1 vector)
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu),
%                 2=superbasic (between bounds), 3=basic (between bounds).
% xs        Initial x vector (nx1), optionally including m slacks at the end.
%           If Warm start, full n+m vector xs must be supplied.
% nS        # of superbasics. Used if a Warm Start, otherwise set to 0.
% SpecsFile Name of the SPECS input parameter file, TOMLAB /SNOPT guide.
% PrintFile Name of the Print file. Name includes the path, maximal number 
%           of characters = 500
% SummFile  Name of the Summary file. Name includes the path, maximal number
%           of characters = 500
% ObjAdd    constant added to the objective for printing purposes, typically 0.
% moremem   Add extra memory for the sparse LU, might speed up the optimization.
%           1E6 is 10MB of memory. If empty, set to 0.
% ProbName  Name of the problem. <=100 characters are used in the MEX interface.
%           In the SQOPT solver the first 8 characters are used in the printed
%           solution. Blank is OK.
% Prob      Sending the Prob structure is optional, only of use if sending
%           H as a function string, see input H.
%
% DESCRIPTION of optPar vector:
%
% Use missing value (-999 or less), when no change of parameter setting is
% wanted. The default value will then be used by SQOPT, unless the value is
% altered in the SPECS file (input SpecsFile). Refer to sqoptTL.m for 
% information about optPar settings.
%
% ------------------------------------------------------------------------
%
% OUTPUT: 
% xs        Solution vector (n+m by 1) with n decision variable values 
%           together with the m slack variables.
% hs        Basis status of variables + constraints (n+m x 1 vector).
%           State of variables: 0=nonbasic (on bl), 1=nonbasic (on bu),
%                 2=superbasic (between bounds), 3=basic (between bounds).
% pi        Lagrangian multipliers (dual solution vector) (m x 1 vector).
% rc        Reduced costs, a n+m vector. If nInf=0, last m = pi.
% Inform    Result of SQOPT run.
%           0 = Optimal solution found.
% nS        # of superbasics.
% nInf      Number of infeasibilities.
% sInf      Sum of infeasibilities.
% Obj       Objective function value at optimum.
% iwCount   Number of QP iterations in iwCount(1), number of Hx products.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written July 16, 2000.  Last modified Dec 24, 2005.
%# mex

function ... 
[xs, hs, pi, rc, Inform, nS, nInf, sInf, Obj, iwCount] = ...
     sqopt( A, bl, bu, H, c, hElast, Prob, iObj, optPar, ...
            Warm, hs, xs,  nS, ...
            SpecsFile, PrintFile, SummFile, ...
            PriLev, ObjAdd, moremem, ProbName )

help sqopt;
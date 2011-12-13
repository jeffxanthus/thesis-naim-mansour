% function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,Prob)
%
% CPLEX MIP User cut callback.
%
% The User Cut Callback is enabled by setting callback(15)=1 in the call to
% cplex.m, or Prob.MIP.callback(15)=1 if using tomRun('cplex',...)
%
% This callback is called by CPLEX during MIP branch & cut for every 
% node that has an LP optimal solution with objective value below the 
% cutoff and is integer infeasible. CPLEX also calls the callback when 
% comparing an integer feasible solution, including one provided by a 
% MIP start before any nodes exist, against lazy constraints.
%
% The callback routine can add globally valid cuts to the LP subproblem. A cut
% is a constraint of the following form:
%
%    c1*x(1) + c2*x(2) + ... + cn*x(n) <?>  rhs
% 
% where <?> is exactly one of the relations <=, >= or = and rhs is a scalar
% right hand side limit. 
%
% Calling syntax:
%
% function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,Prob)
%  
% cpxcb_USERCUT is called by the solver with four arguments:
%
%  x     - The new integer solution
%  f     - The objective value at x
%  xstat - Status information for each variable.  
%  Prob  - The TOMLAB problem structure
%
% cpxcb_USERCUT should return four values: 
%
% ret   - Should return one of the following scalar values:
%
%  0    - Continue optimization
%  1    - Terminate optimization
%         Any other return value will be interpreted as 0.
%
% cut   - A sparse or dense matrix with cuts to add to the current LP
%         subproblem. This matrix must be a m*n sparse or dense double matrix,
%         where m is the number of cuts to add and n is the number of variables
%         (equal to length(x)). If no cuts are to be added, return [] (empty).
%
% sense - If the 'cut' return argument is nonempty, this array should be a m*1 
%         character array indicating the sense of each cut. Allowed values are 
%         L, <, G, >, E, =, for 'less than', 'greater than' and 'equal' respectively. 
%         If more than one cut is specified, simply stack the values, e.g.
%         'LEE'.
%
% rhs   - If the 'cut' return argument is nonempty, this array should be a m*1
%         array of right hand side values for each cut. Cuts are single-sided
%         constraints, thus only one rhs value exists for each cut.
%
% If modifying this file, it is recommended to make a copy of it which
% is placed before the original file in the MATLAB path.
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2008 by Tomlab Optimization Inc., $Release: 11.2.0$
% Written Aug 1, 2007.  Last modified Dec 10, 2008.

function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,Prob)

% NOTE: Important that empty outputs are returned if no cuts are added. 
cut = []; sense = []; rhs = [];

% User code here to optionally set values in cut, sense, rhs. 

% Return value: 0 to continue, 1 to terminate CPLEX. 
ret = 0;

% MODIFICATION LOG:
%
% 070801 ango Wrote function

% function ret = cpxcb_PRESOLVE(cpxCBInfo)
%
% CPLEX Presolve callback
%
% Called at regular intervals from TOMLAB /CPLEX during presolve.
%
% This callback is enabled by setting callback(6)=1 in the call to
% cplex.m, or Prob.MIP.callback(6)=1 if using tomRun('cplex',...)
%
% cpxcb_PRESOLVE is called with one argument, the cpxCBInfo
% progress information vector. 
%
% Contents of cpxCBInfo variable:
%
%  i  cpxCBInfo(i)    - Value
%  -------------------------------------------------------------
%  1  PRESOLVE_ROWSGONE - number of rows eliminated          
%  2  PRESOLVE_COLSGONE - number of columns eliminated       
%  3  PRESOLVE_AGGSUBST - number of aggregator substitutions 
%  4  PRESOLVE_COEFFS   - number of modified coefficients    
%
% By returning a nonzero value from cpxcb_PRESOLVE, the user can
% terminate the optimization. 
%
% If modifying this file, it is recommended to make a copy of it which
% is placed before the original file in the MATLAB path.
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 11.0.0$
% Written Sept 22, 2002   Last modified Feb 2, 2006

function ret = cpxcb_PRESOLVE(cpxCBInfo)

%% INSERT USER CODE HERE.
%
% To terminate optimization, return a nonzero value in 'ret'

ret = 0;

% function ret = cpxcb_MIPPROBE(cpxCBInfo)
%
% CPLEX MIP Probe and Clique Merging callback
%
% Called from TOMLAB /CPLEX during MIP Probing and Clique Merging
%
% This callback is enabled by setting callback(8)=1 in the call to
% cplex.m, or Prob.MIP.callback(8)=1 if using tomRun('cplex',...)
%
% cpxcb_MIPPROBE is called with one argument, the cpxCBInfo progress
% information vector.
%
% Contents of cpxCBInfo variable:
%
%  i  cpxCBInfo(i)    - Value
% -------------------------------------------------------------
%  1 BEST_INTEGER     - obj. value of best integer solution						       
%  2 BEST_REMAINING   - obj. value of best remaining node 						       
%  3 NODE_COUNT       - total number of nodes solved  						   	       
%  4 NODES_LEFT       - number of remaining nodes 							       
%  5 MIP_ITERATIONS   - total number of MIP iterations 						       
%  6 MIP_FEAS         - returns 1 if feasible solution exists; otherwise, 0 				       
%  7 CUTOFF           - updated cutoff value 								       
%  8 CLIQUE_COUNT     - number of clique cuts added 							       
%  9 COVER_COUNT      - number of cover cuts added 							       
% 10 DISJCUT_COUNT    - number of disjunctive cuts added 						       
% 11 FLOWCOVER_COUNT  - number of flow cover cuts added 						       
% 12 FLOWPATH_COUNT   - number of flow path cuts added 						       
% 13 FRACCUT_COUNT    - number of Gomory fractional cuts added 					       
% 14 GUBCOVER_COUNT   - number of GUB cover cuts added 						       
% 15 IMPLBD_COUNT     - number of implied bound cuts added 						       
% 16 MIRCUT_COUNT     - number of mixed integer rounding cuts added 					       
% 17 PROBE_PHASE      - current phase of probing (0-3) 						       
% 18 PROBE_PROGRESS   - fraction of probing phase completed (0.0-1.0) 				   	       
% 19 FRACCUT_PROGRESS - fraction of Gomory cut generation for the pass completed (0.0 - 1.0)  	   	       
% 20 DISJCUT_PROGRESS - fraction of disjunctive cut generation for the pass completed (0.0 - 1.0)  	       
% 21 FLOWMIR_PROGRESS - fraction of flow cover and MIR cut generation for the pass completed (0.0 - 1.0)      
% 22 MY_THREAD_NUM    - identifier of the parallel thread making this call (always 0)
% 23 USER_THREADS     - total number of parallel threads currently running (always 1)
%
% By returning a nonzero value from cpxcb_MIPPROBE, the user can
% terminate the optimization. 
%
% If modifying this file, it is recommended to make a copy of it which
% is placed before the original file in the MATLAB path.
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 11.0.0$
% Written Sept 22, 2002   Last modified Feb 2, 2006

function ret = cpxcb_MIPPROBE(cpxCBInfo)

%% ADD USER CODE HERE.
%
% To terminate optimization return a nonzero value in 'ret'

ret = 0;

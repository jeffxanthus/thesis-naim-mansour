% Initialization of structure Result
%
% function Result = ResultDef(Prob)
%
% INPUT:
%  Prob     Structure for the problem
%
% OUTPUT:
%  Result   Structure for the results solving problem Prob.P

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written May 14, 1998.   Last modified Jun 6, 2008.

function Result = ResultDef(Prob)

global probType

solvType=Prob.solvType;

Result = struct('Name',Prob.Name, 'P',Prob.P, 'probType',probType, 'Solver','', ...
  'SolverAlgorithm','', 'solvType',solvType, 'ExitFlag',[], 'ExitText',[], ...
  'Inform',[], 'CPUtime','', 'REALtime','', ...
  'Iter',[], 'MinorIter',double(0), 'maxTri',[], ...
  'FuncEv',double(0), 'GradEv',double(0), 'HessEv',double(0), ...
  'ConstrEv',double(0), 'ConJacEv',double(0), 'ConHessEv',double(0), ...
  'ResEv',double(0), 'JacEv',double(0), ...
  'x_k',[], 'f_k',[], 'g_k',[], 'B_k',[], 'H_k',[],  ...
  'y_k',[], 'v_k',[], 'r_k',[], 'J_k',[], 'Ax',[], 'c_k',[], 'cJac',[], ...
  'd2L',[], ...
  'x_0',[], 'f_0',[], 'c_0',[], 'Ax0',[],  ...
  'xState',[], 'bState',[], 'cState',[], 'p_dx',[],  'alphaV',[], ...
  'x_min',[], 'x_max',[], 'LS', [], 'F_X',[], 'SepLS',[], ...
  'QP',[], 'SOL',[]);

% MODIFICATION LOG:
%
% 980919  hkh  New fields due to change from optPar to structure optParam
% 980922  hkh  New field H_k, either Quasi-Newton matrix or Hessian
% 981011  hkh  Deleted field optParam. Prob.optParam is used instead.
%              Added fields CPUtime, REALtime, Nflops
% 981026  hkh  Add field SolverAlgorithm, a text description of the algorithm
% 981028  hkh  Add field f_0, initial function value
% 981108  hkh  Add field SepLS, for separable NLLS
% 981111  hkh  One field for Hessian, H_k, one for Quasi-Newton matrix, B_k
% 000708  hkh  Add MinorIter, for iterations in sub problem
% 000721  hkh  Add SOL field
% 000916  hkh  Add ExitText field
% 001009  hkh  Removed Nflops field, flops removed in 6.0
% 011031  hkh  Added maxTri field for global optimization with rectangles
% 020811  hkh  Added LS field for output statistics
% 030309  hkh  Added field Ax for evaluation of linear equations A*x
% 040101  hkh  Use double(0), also for MinorIter
% 040109  hkh  Rearranging fields, add fields c_0 and Ax0
% 040407  hkh  Adding counter fields ConJacEv and ConHessEv
% 050223  frhe Adding d2L_k to Result structure
% 080606  med  isfield checks removed
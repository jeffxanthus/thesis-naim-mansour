% odeSolverList returns a list of all ordinary differential equations solvers
%
% function SolvList =  odeSolverList(Silent)
%
% INPUT:
%
% Silent       If Silent == 1, totally silent
%              If Silent == 0 (default), an information text is displayed
%
% OUTPUT:
% SolvList     A list of all solvers for ordinary differential equations

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function [SolvList] = odeSolverList(Silent)

solveType = 1;
nSolvTypes = 1;

if nargin < 2
   Silent = 0;
   if nargin < 1
      solvType = [];
   end
end

SolvList=str2mat('lsode','rksuite','modfit');
RecSolver='lsode';

if Silent, return; end

% Print recommended choices of solvers

fprintf('\nTomlab recommended choice for ');
fprintf('Ordinary differential equations (ODE)');
fprintf('\n\n');
fprintf(RecSolver);
fprintf('\n\n');

% MODIFICATION LOG
%
% 050412   bkh   Written from copy of SolverList
% 050705   med   Help updated
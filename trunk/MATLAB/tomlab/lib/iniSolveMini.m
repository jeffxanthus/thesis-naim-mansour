% init routine. Called before solving an optimization problem
%
% Initializes global variables
%
% function Prob = iniSolveMini(Prob,solvType);
%
% solvType   Solver type

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Jun 6, 2008.  Last modified Feb 5, 2009.

function Prob = iniSolveMini(Prob)

if nargin < 1
    error('Prob not supplied');
end

if Prob.smallA
   % Find small elements in A and remove them
   v         = max(max(abs(Prob.A)));
   [i,j,s]   = find(sparse(Prob.A));
   ix        = abs(s) < eps*v;
   A0        = sum(ix);
   if A0 > 0
      [m,n]  = size(Prob.A);
      if Prob.Warning | Prob.optParam.IterPrint | Prob.PriLevOpt > 0 
         fprintf('iniSolve: Delete %d small elements in A matrix\n',A0);
      end
      ix1    = find(ix == 0);
      % Create new A matrix
      Prob.A = sparse(i(ix1),j(ix1),s(ix1),m,n);
   end
end

Prob.TIME0 = cputime;
Prob.TIME1 = clock;

% Communication qp_c/dc - reset for CPLEX
global QP_Qx QP_Qxx
QP_Qx = []; QP_Qxx = [];

% MODIFICATION LOG
%
% 080606  med  Wrote file
% 080607  hkh  Add PreSolve with comments, remove global variables
% 090205  med  Global reset for miqq problems
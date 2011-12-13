% odeRun - General driver routine for a TOMLAB ODE problem
%
% INPUT
%
%  Solver     The name of the ODE solver used to solve the ODE system
%  Prob       TOMLAB problem structure
%  PriLev     Print level when displaying the result
%
% OUTPUT
%
%  Result     Tomlab result structure

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Apr 21, 2005.  Last modified Aug 02, 2005.

function Result = odeRun(Solver, Prob, PriLev)

global n_ode_f  n_ode_J  % Count of function and Jacobian evaluations

if nargin < 3
   PriLev = [];
   if nargin < 2
      error([' odeRun must have at least two parameters, Solver and Prob' ...
           ' structure']);
   end
end

if ~isstruct(Prob)
   error('Input parameter Prob must beof type struct');
end

if isempty(Solver)
   Solver = 'lsode';
elseif ~ischar(Solver)
   Solver = 'lsode';
end

% Check lower level ODE parameters in Prob.ODE
Prob = odeProbCheck(Prob);
 
switch lower(Solver)
 case 'lsode'
   Result            = lsodeTL(Prob);
 case {'dopri5', 'radau5', 'ind-dir'}
   Result            = modfitodeTL(Prob, Solver);
 case {'ode23' ,'ode113' ,'ode15i' ,'ode15s' ,'ode23s' ,'ode23t' ...
      ,'ode23tb' ,'ode45'}
   Result            = odeML(Prob, Solver);
 case 'rksuite'
   Result            = rksuiteTL(Prob);
 otherwise
   fprintf('Solver %s', Solver);
   fprintf(' NOT found\n');
   error('Illegal solver algorithm!');
end
Result.ODE.FuncEv = n_ode_f;
Result.ODE.JacEv  = n_ode_J;

% MODIFICATION LOG:
%
% 050413  joho  Written
% 050415  joho  Changed global variables n_f, n_J to n_ode_f, n_ode_J
% 050417  hkh   Cleaned up, changed funcev to FuncEv
% 050429  bkh   Added matlab and modfit solvers
% 050502  hkh   Minor changes.
% 050502  hkh   Add call to odeProbCheck, ODE solution specific fields
% 050705  med   Help updated
% 050801  med   isstr replaced by ischar
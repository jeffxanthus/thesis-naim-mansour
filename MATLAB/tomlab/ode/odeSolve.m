% function Result = odeSolve(odeSolver,Prob)
%
% odeSolve solves an ODE sub problem without any check on the Prob structure
%
% It is intended for use by other solvers when solving an ODE subproblem
%
% INPUT PARAMETERS
% odeSolver  Name of the ODE solver
% Prob       Input structure, feeded to the solver
%
% OUTPUT PARAMETERS
% Result     Output result structure, feeded back to the caller.

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Result = odeSolve(odeSolver,Prob)

if nargin < 2
   error('odeSolve needs input structure Prob');
end

if isempty(Prob.ODE.J)
   jac = [];
else
   jac = 'ode_J';
end

switch lower(odeSolver)
   case {'lsode'}
      [Result.ODE.y, Result.Inform] = lsode(Prob.ODE.Y0, Prob.ODE.tInit, ...
      Prob.ODE.tWant, 'ode_f', Prob.ODE.J, Prob.ODE.LSODE, Prob.PriLevOpt, Prob);
   case {'rksuite'}
      [Result.ODE.y, Result.Inform] = rksuite(Prob.ODE.Y0, Prob.ODE.tInit, ...
      Prob.ODE.tStop, Prob.ODE.tWant, Prob.ODE.InitStep, Prob.ODE.f, ...
      Prob.PriLevOpt, Prob.ODE.RKSUITE, Prob);
   case {'dopri5','radau5','ind-dir'}
       % Temporary solution
       Result = modfitodeTL(Prob, odeSolver);
   case {'ode45'} 
      [t, Result.ODE.y] = ode45('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode23'} 
      [t, Result.ODE.y] = ode23('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode113'} 
      [t, Result.ODE.y] = ode113('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode15i'} 
      [t, Result.ODE.y] = ode15i('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode15s'} 
      [t, Result.ODE.y] = ode15s('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode23s'} 
      [t, Result.ODE.y] = ode23s('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode23t'} 
      [t, Result.ODE.y] = ode23t('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   case {'ode23tb'} 
      [t, Result.ODE.y] = ode23tb('odeML_f', Prob.ODE.tWant, Prob.ODE.Y0, ...
                                 Prob.ODE.ML, Prob);
   %case {'ode23' ,'ode113' ,'ode15i' ,'ode15s' ,'ode23s' ,'ode23t' ...
   %     ,'ode23tb' ,'ode45'}
   %   [t, Result.ODE.y] = feval(lower(odeSolver),'odeML_f', Prob.ODE.tWant, ...
   %   Prob.ODE.Y0, Prob.ODE.ML, Prob);
   otherwise
      Result = odeRun(odeSolver,Prob);
end

% MODIFICATION LOG:
%
% 050421  bkh  Created from copy of tomSolve
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to
%              InitStep, tInit and tStop respectively
% 050429  bkh  Added matlab and modfit odesolvers
% 050502  hkh  Minor changes in comments. Avoid feval for Matlab ODEs
% 050705  med  Help updated
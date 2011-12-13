% odeML - Interface for solving an ODE problem set in the TOMLAB standard
%         format using MATLABs internal ODE solvers.
%
% INPUT:
%
%  Prob       Tomlab problem structure
%  Solver     The name of the MATLAB ODE solver used to solve the ODE system 
%
% -------------------------------------------------------------------------
% Fields used in input structure Prob 
% -------------------------------------------------------------------------
% PriLev    Print level.
%
% -------------------------------------------------------------------------
% Fields used in Prob.ODE:
% -------------------------------------------------------------------------
% f         Name of the function of the ODE system, f(t,y).
% J         Name of the Jacobian of the ODE system, J(t,y).
% Y0        The initial values of the ODEs. nEq = length(Y0) is the number
%           of ordinary differential equations for the problem
% tInit     Initial values of the independent variables
% tWant     Values of the independent variables to return function values for
% tStop     The solver integrates from tInit to tStop if tWant is not given
% relTol    Relative tolerance parameter. Given either as a
%           scalar or a vector of length nEq. The minimum 
%           relative tolerance is 1e-16 and the maximum is 1e-1
% absTol    Absolute tolerance parameter. Given either as a 
%           scalar or a vector of length nEq
% InitStep  The step size to be attempted on the first step.
%           Default value is determined by the solver
% ML        Structure containing the following extra options:
%
%           NormControl
%           Refine
%           Stats
%           JPattern
%           MaxStep
%           BDF
%           MaxOrder
%
%           See help odeset for information on these options
%
% -------------------------------------------------------------------------
%
% OUTPUT:
%
%  Result     Tomlab result structure
%
% ODE.t    Time vector t where ODE intgreation values y are reported
% ODE.y    Matrix length(t) x nEq of integration values y at the time points t

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Result = odeML(Prob, Solver)

if nargin < 2
   Solver = 'ode45';
   if nargin < 1 
      error('odeML needs the Prob structure as input');
   end
end

odeFun = 'odeML_f';
tSpan  = DefPar(Prob.ODE, 'tWant', []);
if isempty(tSpan)
   tStop = DefPar(Prob.ODE, 'tStop', []);
   if isempty(tStop)
      error('Neither tWant nor tStop set!');
   else
      tSpan = [tInit tStop];
   end
end
Y0     = Prob.ODE.Y0;

% Set fields in struct Options corresponding to fields in Prob.ODE
% and fields in Prob.ODE.ML
Prob.ODE.ML    = DefPar(Prob.ODE, 'ML', []);
Options        = Prob.ODE.ML;
Options.RelTol = DefPar(Prob.ODE, 'relTol', []);
Options.AbsTol = DefPar(Prob.ODE, 'absTol', []);
if isfield(Prob.ODE, 'InitialStep')
   Options.InitialStep = Prob.ODE.InitStep;
end
if isempty(Prob.ODE.J)
   % no jacobian
else
   Options.Jacobian = 'odeML_J';
end

switch Solver
  case 'ode23'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode113'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode15i'
    error('Solver ode15i is not supported yet!');
  case 'ode15s'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode23s'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode23t'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode23tb'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  case 'ode45'
    [Result.ODE.t,Result.ODE.y] = ode45(odeFun, tSpan, Y0, Options, Prob);
  otherwise
    disp('Illegal choice of odesolver!');
end

% MODIFICATION LOG
%
% 050418  bkh  Created modification log to this new file
% 050418  bkh  Supports ode45
% 050419  bkh  Supports all solvers except ode15i
% 050419  bkh  Added ML specific options
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050502  hkh  Minor comment changes, clean up
% 050705  med  Help updated
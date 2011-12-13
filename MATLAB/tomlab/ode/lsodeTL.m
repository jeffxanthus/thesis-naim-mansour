% TOMLAB LSODE Ordinary Differential Equations Solver
%
% function Result = lsodeTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format.
%
% -------------------------------------------------------------------------
% Fields used in input structure Prob
% -------------------------------------------------------------------------
% PriLev    Print level.
%
% -------------------------------------------------------------------------
% Fields used in Prob.ODE:
% -------------------------------------------------------------------------
%
% f         Name of the function of the ODE system, f(t,y).
% J         Name of the jacobian of the ODE system, jac(t,y).
% Y0        The initial values of the ODEs. The length of y is the number
%           of equations in the system.
% tInit     Initial values of the independent variables.
% tWant     Values of the independent variables to return function values for.
% tStop     Critical value of t which the solver is not to overshoot.
%           Only required if LSODE.Task is 4 or 5.
% relTol    Relative tolerance parameter. Given either as a
%           scalar or a vector of length neq. The minimum
%           relative tolerance is 1e-16 and the maximum is 1e-1.
% absTol    Absolute tolerance parameter. Given either as a
%           scalar or a vector of length neq.
% InitStep  The step size to be attempted on the first step.
%           Default value is determined by the solver.
% LSODE     Structure containing the following extra options:
%           Task      Specifies the task
%                      1 - Normal computation of y(t) by overshooting and
%                          interpolation.
%                      2 - Take one step only and return
%                      3 - Stop at the first internal mesh point at or
%                          beyond t and return.
%                      4 - Normal computation of output values of y(t) at
%                          t but without overshooting t = tCrit.
%                          tCrit must be given in the LSODE structure.
%                      5 - Take one step without passing tCrit and return.
%                          tCrit must be given in the LSODE structure.
%           Method    Specifies the method
%                      1 - The implicit Adams method.
%                      2 - Method based on backward differentiation
%                          formulas (bdf-s).
%           JacType   Set what kind of jacobian should be used.
%                     Six different options are availible.
%                      0 - No internally generated jacobian.
%                      1 - Internally generated diagonal jacobian.
%                      2 - Internally generated banded jacobian.
%                      3 - Internally generated full jacobian.
%                      4 - User supplied banded jacobian
%                      5 - User supplied full jacobian
%           JacLHBW   The lower half bandwith of a user supplied banded
%                     jacobian.
%           JacUHBW   The upper half bandwith of a user supplied banded
%                     jacobian.
%                     The elements within the band of the jacobian should
%                     be loaded into the returned matrix in a columnwise
%                     manner, with the diagonal lines loaded into rows.
%                     Thus df(i)/dy(i) is loaded into position
%                     (i-j+JacUHBW+1, j).
%           MaxOrder  The maximum order to be allowed. The default value is
%                     12 if Method = 1, and 5 if Method = 2.
%                     The default values are also the maximum values.
%           MaxSteps  The maximum number of allowed steps for each
%                     computation at one value of the variable t.
%                     The default value is 500.
%           hMax      The maximum absolute step size allowed.
%                     The default value is infinite.
%           hMin      The minimum absolute step size allowed.
%                     The default value is 0. This lower bound is not
%                     enforced on the final step before reaching tCrit
%                     when Task = 4 or 5.)
%           MaxhWarn  Maximum number of messages printed (per problem)
%                     warning that t + h = t on a step (h = step size).
%                     this must be positive to result in a non-default
%                     value. The default value is 10.
%
% -------------------------------------------------------------------------
%
% OUTPUT:
%
% Result     Structure with results (see ResultDef.m):
%
% The following outputs are set in the Result.ODE sub field:
%
% y          Solution at t.
% t          Value of independent value in solution y.
%
% The following outputs are set in the Result.ODE.LSODE sub field:
%
% Inform     2 - lsode was successful.
%            1 - Nothing to do, tWant = tInit.
%           -1 - Excess work done on this call.
%           -2 - Excess accuracy requested (tolerances too small).
%           -3 - Illegal input detected (see printed message).
%           -4 - Repeated error test failures (check all inputs).
%           -5 - Repeated convergence failures (perhaps bad jacobian
%                supplied or wrong choice of tolerances).
%           -6   Error weight became zero during problem. (solution
%                component i vanished, and absTol or absTol(i) = 0.)
% tReach    The independent variable. Usually the same as the last value
%           in the vector t, but in case of error, tReach is the farthest
%           point reached.
% hUsed     A vector with the sizes of each successful step in t
% TolScalF  A tolerance scale factor, greater than 1.0,
%           computed when a request for too much accuracy was
%           detected (Inform = -3 if detected at the start of
%           the problem, Inform = -2 otherwise). If relTol and
%           absTol are uniformly scaled up by a factor of
%           TolScalF for the next call, then the solver is
%           deemed likely to succeed. (the user may also ignore
%           TolScalF and alter the tolerance parameters in any
%           other way appropriate.)
% nSteps    The number of steps taken for the problem so far.
% nFeval    The number of f evaluations for the problem so far.
% nJeval    The number of jacobian evaluations (and of matrix
%           lu decompositions) for the problem so far.
%
% -----------------------------------------------------------------------
%
% For a problem description, see odeAssign.m
%
% -------------------------------------------------------------------------

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Result = lsodeTL(Prob)

if nargin < 1 
    error('lsodeTL needs the Prob structure as input');
end

Y0               = DefPar(Prob.ODE      , 'Y0'      , []);
tInit            = DefPar(Prob.ODE      , 'tInit'   , []);
tWant            = DefPar(Prob.ODE      , 'tWant'   , []);
f                = 'ode_f';
jac              = 'ode_J';

Options          = DefPar(Prob.ODE      , 'LSODE'   , []);
Options.relTol   = DefPar(Prob.ODE      , 'relTol'  , []);
Options.absTol   = DefPar(Prob.ODE      , 'absTol'  , []);
Options.Task     = DefPar(Options       , 'Task'    , 1);
if Options.Task > 3
   Options.tStop = DefPar(Prob.ODE      , 'tStop'   , []);
end
Options.InitStep = DefPar(Prob.ODE      , 'InitStep', []);
PriLev           = DefPar(Prob.PriLevOpt, 'PriLev'  , []);

[yOut, Inform, tOut, tReach, hUsed, TolScalF, nSteps, nFeval, nJeval] = ...
lsode(Y0, tInit, tWant, f, jac, Options, PriLev, Prob);

if min(Inform) < 0
   Result.Inform = min(Inform);
   Result.ExitFlag = 0;
else
   Result.Inform = max(Inform);
   Result.ExitFlag = 1;
end
switch Result.Inform
    case -6
        Result.ExitText = 'Error weight became zero during problem. (solution component i vanished, and AbsTol or AbsTol(i) = 0.) ';
    case -5
        Result.ExitText = 'Repeated convergence failures (perhaps bad jacobian supplied or wrong choice of tolerances).';
    case -4
        Result.ExitText = 'Repeated error test failures (check all inputs).';
    case -3
        Result.ExitText = 'Illegal input detected (see printed message).';
    case -2
        Result.ExitText = 'Excess accuracy requested (tolerances too small).';
    case -1
        MaxSteps = 500;
        if isfield(Prob.ODE.LSODE, 'MaxSteps')
            if ~isempty(Prob.ODE.LSODE.MaxSteps)
               MaxSteps = Prob.ODE.LSODE.MaxSteps;
            end
        end                    
        Result.ExitText = ['Excess work done. The number of steps exceeds MaxSteps = ', num2str(MaxSteps), ...
                           '. Change the parameter MaxSteps and/or reinitialize at the last successful t value, t = ', num2str(t)];
    case  1
        Result.ExitText = 'Nothing to be done, tWant = tInit';
    case  2
        Result.ExitText = 'Complete success!';
end

Result.ODE.LSODE.Inform   = Inform;
Result.ODE.y              = yOut;
Result.ODE.t              = tOut;
Result.ODE.LSODE.tReach   = tReach;
Result.ODE.LSODE.hUsed    = hUsed;
Result.ODE.LSODE.TolScalF = TolScalF;
Result.ODE.LSODE.nSteps   = nSteps;
Result.ODE.LSODE.nFeval   = nFeval;
Result.ODE.LSODE.nJeval   = nJeval;

% MODIFICATION LOG
%
% 050418  bkh  Created modification log to this file
% 050418  bkh  Added help, removed commented lines
% 050418  bkh  Change name of yStart to Y0 in Prob.ODE
% 050419  bkh  Help corrected
% 050421  bkh  Changed output order from [Inform, y,...] to [y, Inform,..]
% 050421  bkh  Corrected and simplified setting of struct Options
% 050421  bkh  Set Task = 1 as default, don't set tCrit unless Task = 4 or 5
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050705  med  Help updated
% LSODE Matlab Solver
% -----------------------------------------------------------------------
%
%   lsode solves the following initial value problem for stiff or nonstiff
%   systems of first order ordinary differential equations,
%
%          dy/dt = f(t,y), or, in component form,
%          dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
%
% -----------------------------------------------------------------------
%
% function [yOut, Inform, tOut, t, hUsed, TolScalF, nSteps, nFeval, nJeval]
%          = lsode(Y0, tInit, tWant, f, jac, neq, options, PriLev, Prob)
%
% -----------------------------------------------------------------------
% INPUT:  Must give at least 4 first parameters, the rest are optional.
%
% Y0        The initial values of the ODEs. The length of Y0 is the number.
%           of equations in the system.
% tInit     Initial values of the independent variables.
% tWant     Values of the independent variables to return function values for.
% f         Name of the function of the ODE system, f(t,y).
% jac       Name of the jacobian of the ODE system, jac(t,y).
% options   Structure containing extra options.
% PriLev    Print level, either 0 or 1.
% Prob      Problem structure in TOMLAB format.
%
% -----------------------------------------------------------------------
% Fields used in options:
% -----------------------------------------------------------------------
%
%           relTol    Relative tolerance parameter. Given either as a
%                     scalar or a vector of length neq. The minimum 
%                     relative tolerance is 1e-16 and the maximum is 1e-1.
%           absTol    Absolute tolerance parameter. Given either as a 
%                     scalar or a vector of length neq.
%           Task      Specifies the task
%                      1 - Normal computation of y(t) by overshooting and
%                          interpolation.
%                      2 - Take one step only and return
%                      3 - Stop at the first internal mesh point at or
%                          beyond tWant and return.
%                      4 - Normal computation of output values of y(t) at
%                          t but without overshooting t = tStop.
%                          tStop must be given in the options structure.
%                      5 - Take one step without passing tStop and return.
%                          tStop must be given in the options structure.
%           tStop     Critical value of t which the solver is not to 
%                     overshoot.  Required if Task is 4 or 5.
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
%           InitStep  The step size to be attempted on the first step.
%                     Default value is determined by the solver.
%           hMax      The maximum absolute step size allowed.
%                     The default value is infinite.
%           hMin      The minimum absolute step size allowed.
%                     The default value is 0. This lower bound is not
%                     enforced on the final step before reaching tStop
%                     when Task = 4 or 5.)
%           MaxhWarn  Maximum number of messages printed (per problem)
%                     warning that t + h = t on a step (h = step size).
%                     this must be positive to result in a non-default
%                     value. The default value is 10.
%
% -----------------------------------------------------------------------
% OUTPUT: 
%
% yOut      Solution at t=tOut.
% Inform     2 - lsode was successful.
%            1 - Nothing to do, t = t0.
%           -1 - Excess work done on this call.
%           -2 - Excess accuracy requested (tolerances too small).
%           -3 - Illegal input detected (see printed message).
%           -4 - Repeated error test failures (check all inputs).
%           -5 - Repeated convergence failures (perhaps bad jacobian
%                supplied or wrong choice of tolerances).
%           -6   Error weight became zero during problem. (solution
%                component i vanished, and atol or atol(i) = 0.)
% tOut      Value of independent value in solution yOut.
% t         The independent variable. Usually the same as the last value
%           in the vector tOut, but in case of error, t is the farthest 
%           point reached.
% hUsed     A vector with the sizes of each successful step in t
% TolScalF  A tolerance scale factor, greater than 1.0, 
%           computed when a request for too much accuracy was
%           detected (Inform = -3 if detected at the start of
%           the problem, Inform = -2 otherwise). If rtol and 
%           atol are uniformly scaled up by a factor of 
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

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.
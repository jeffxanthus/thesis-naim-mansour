% RKSUITE ODE TOMLAB Solver
% -----------------------------------------------------------------------
%
%   rksuite is a suite of codes that solves the following initial
%   value problem for a first order system of ordinary differential
%   equations:
%
%      y' = f(t,y),    y(t0) = y0
%
%   where
%
%   y is a nx1 vector of solution components
%   t is the independent (time) variable
%
% -----------------------------------------------------------------------
%
% function ...
% [y, t, yp, yMax, InfLocal, InfGlobal] =
%     rksuite(Y0, tInit, tStop, tWant, InitStep, f, PriLev, RKSUITE, ...
%     Prob);
%
% -----------------------------------------------------------------------
%
% INPUT:  Must give at least 3 first parameters
%
% Y0        Vector of initial values to the ODE system. Sets the
%           dimension of the ODE system.
% tInit     Initial value of the time variable.
% tStop     Integration proceeds from tInit to tStop and rksuite
%           cannot integrate past tStop.
% tWant     Vector with t-values in which rksuite will compute
%           solutions. tWant may be left empty when complicated task
%           is being used.
% InitStep  Initial step length, sets the scale of the problem.
% f         File with the ODE system on the form y' = f(t,y,Prob).
% PriLev    Print Level for rksuite.
% RKSUITE   Matlab structure with rksuite specific settings.
% Prob      Tomlab problem structure.
%
% MEMBERS of input structure RKSUITE:
%
% METHOD    Defines which Runge-Kutta formula pair that rksuite will use.
% TASK      Usual or Complicated task.
% ERRASS    Do or do not try to assess the true error
%
% OUTPUT:
% y         Matrix with solutions corresponding to the values in t.
% t         The values where a solution actually could be found.
% yp        Approximations to the first derivative of the solutions
%           in t.
% yMax      yMax(L) is the largest value of abs(y(L)) at any step
%           of the integration (so far).
% InfLoc    Vector of the rksuite inform value for each value in t.
% InfGlob   Some kind of worst value taken over all values in
%           infLoc.

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

 function ...
 [y, t, yp, yMax, InfLocal, InfGlobal] = ...
     rksuite(Y0, tInit, tStop, tWant, InitStep, f, PriLev, RKSUITE, ...
	     Prob)
 
 help rksuite;
 
% MODIFICATION LOG:
%
% 050419  joho Written
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to
%              InitStep, tInit and tStop respectively
% 050705  med  Help updated
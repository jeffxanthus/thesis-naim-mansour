% function f = odeML_f(t, y, flag, Prob)
%
% TOMLAB gateway routine for the computation of the function values f(y,t)
% using MATLABs internal ODE solvers
%
% odeML_f calls the routine Prob.ODE.f as
%           f=feval(Prob.ODE.f, t, y, Prob)
%
% The global counter variable n_ode_f is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function f = odeML_f(t, y, flag, Prob)

global n_ode_f 

n_ode_f = n_ode_f+1;

switch flag
 case ''
   f = feval(Prob.ODE.f, t, y, Prob);
 case 'init'
   f = Prob.ODE.Y0;
 otherwise
   error(['Unknown flag ' flag]);
end

% MODIFICATION LOG:
%
% 050418 bkh  Written from copy of ode_f.m
% 050418 hkh  Corrected comments, Prob.ODE.f is used, not Prob.FUNCS.f
% 050705 med  Help updated
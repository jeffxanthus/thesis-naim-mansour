% function f = ode_f(t, y, Prob)
%
% TOMLAB gateway routine for the computation of the function values f(y,t)
%
% ode_f calls the routine Prob.FUNCS.f as 
%           f=feval(Prob.FUNCS.f, t, y, Prob)
%
% The global counter variable n_ode_f is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function f = ode_f(t, y, Prob)

global n_ode_f 

n_ode_f = n_ode_f+1;

f = feval(Prob.ODE.f, t, y, Prob);

% MODIFICATION LOG:
%
% 050407 bjjo Written
% 050417 hkh  Cleaned up, corrected minor errors
% 050418 hkh  Revision, simplified
% 050705 med  Help updated
% function J = ode_J(t, y, Prob)
%
% TOMLAB gateway routine for the computation of the Jacobian J(y,t)
%
% ode_J calls the routine Prob.ODE.J as 
%           J=feval(Prob.FUNCS.J, t, y, Prob)
%
% The global counter variable n_ode_J is incremented

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function J = ode_J(t, y, Prob)

global n_ode_J 

n_ode_J = n_ode_J+1;

f = feval(Prob.ODE.J, t, y, Prob);

% MODIFICATION LOG:
%
% 050407 bjjo Written
% 050415 hkh  Cleaned up, corrected minor errors
% 050418 hkh  Revision, simplified
% 050524 med  Changed counter, said _f before
% 050705 med  Help updated
% function cx=od_c(x, Prob)
%
% The od problems are Constrained Nonlinear Least Squares Problems
% where the function is an ODE
%
% od_c defines the nonlinear constraints

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.8.0$
% Written Apr 18, 2005.  Last modified Apr 18, 2005.

function cx=od_c(x, Prob)

P  = Prob.P;
cx = [];

if P==1
    % 'Arporn Reactor'
end

% MODIFICATION LOG:
%
% 050418  hkh  Defined

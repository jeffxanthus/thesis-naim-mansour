% function r = odeFit_r(x, Prob)
%
% TOMLAB callback routine for the residual vector for an ODE parameter
% estimation problem.
%
% INPUT:
% x       Current iterate x, where the residual is vector function r=r(x).
% Prob    TOMLAB problem structure
%
% OUTPUT:
% r       Non-weighted residual vector

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function r = odeFit_r(x, Prob)

Y0Idx = Prob.ODE.Y0Idx;
if isempty(Y0Idx)
   Prob.ODE.X         = x;
else
   % Split ODE initial value parameters from other unknown parameters
   Prob.ODE.Y0(Y0Idx) = x(1:length(Y0Idx));
   Prob.ODE.X         = x(length(Y0Idx)+1:end);
end

R  = odeSolve(Prob.SolverODE, Prob);
Y  = R.ODE.y(:);

% ODE Model - Data
r = Y(Prob.ODE.yIdx) - Prob.LS.y;

% MODIFICATION LOG:
%
% 050413  joho Written
% 050415  joho Bugfixes and change of documentation
% 050418  hkh  Revised, made more efficient handling of indices in ODE.y
% 050418  hkh  Split x into ODE.Y0 and ODE.X
% 050705  med  Help updated
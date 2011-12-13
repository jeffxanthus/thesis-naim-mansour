% odFunc_f.m
%
% function f = odFunc_f(t, y, Prob)
%
% Computes ODE function used in integrating numerically the ODE

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 18, 2005.   Last modified Apr 18, 2005.

function dy = odFunc_f(t, y, Prob)

P = Prob.P;

% Parameters varied are sent using the Tomlab structure, in Prob.ODE.X
x     = Prob.ODE.X;
% Constant parameters are sent using the Tomlab structure, in Prob.ODE.user
param = Prob.ODE.param;
% Constant vector of parameters, Z,  sent in Prob.ODE.Z

if P==1
    % 'Arporn Reactor'
    valG  = param.valG;
    valL1 = param.valL1;
    valL2 = param.valL2;
    valL3 = param.valL3;
    valT  = param.valT;
    H     = param.H;
    R     = param.R;

    rate1 = x(1)*exp(-x(2)/(R*y(6)))*y(3)*y(2);
    rate2 = x(3)*exp(-x(4)/(R*y(6)))*y(4)*y(2);

    dy    = [ valG*(y(1)-y(2));
        valL1*(rate1 + rate2) + valL2*(y(1)-y(2));
        valL1*(rate1);
        valL1*(rate2 - rate1);
        valL3*(rate2);
        valT*(-H)*(rate1 + rate2)];
end

% MODIFICATION LOG:
%
% 050418  hkh  Written, defined Arporn reactor
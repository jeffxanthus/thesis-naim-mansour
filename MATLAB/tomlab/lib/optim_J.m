% optim_J.m
%
% function J = optim_J(x, Prob)
%
% optim_J is used to implement the OPT TB 2.x interface
%
% The Jacobian matrix J is returned, if available in the global variable LS_J
% with the corresponding x value in LS_xJ
%
% optim_J is called from the TOMLAB gateway function nlp_J.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 29, 1999.  Last modified Jan 29, 2003.

function J = optim_J(x, Prob)

global LS_xJ LS_J

% Return the gradient, if computed.

if ~isempty(LS_xJ)
    if all(x==LS_xJ)
        J = LS_J;
    elseif ~isempty(LS_J)
        r = optim_rJ(x, Prob);
        J = LS_J;
    else
        J = [];
    end
else
    J = [];
end

% MODIFICATION LOG:
%
% THIS ROUTINE IS NORMALLY NOT NEEDED, BECAUSE nlp_J detects that
% LS_J is computed with the correct LS_xJ value
% 030129  hkh  Call optim_rJ, not optim_r
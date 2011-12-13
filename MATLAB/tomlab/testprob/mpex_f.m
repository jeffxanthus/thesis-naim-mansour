% function f = mpex_f(x, Prob)
%
% Compute function value to problem Prob.P
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% f      Function value, f(x).

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Jun 22, 2006.  Last modified Jun 3, 2008.

function f = mpex_f(x,Prob)
P = Prob.P;
if P == 1
    % Gauvin and Savard
    f = x(1)^2 + (x(2)-10)^2;
elseif P == 2
    % DF1
    f = (x(1) - 1 - x(2))^2;
elseif P == 3
    % BILIN
    f = 8*x(1) + 4*x(2) - 4*x(3) + 40*x(4) + 4*x(5);
elseif P == 4
    % Bard3m
    f = -x(1)^2 - 3*x(2) + x(4)^2 - 4*x(3);
elseif P == 5
    % RalphMOD
    P = Prob.user.P;
    c = Prob.user.c;
    nd = 4;
    ns = 100;
    N = nd+ns;
    f = 0.5*( ...
        x(1:nd)' * P(1:nd,1:nd) * x(1:nd) + ...
        + 2.0*x(1:nd)'*P(1:nd,nd+1:N)*x(nd+1:N) + ...
        x(nd+1:N)'*P(nd+1:N,nd+1:N)*x(nd+1:N) ) + ...
        c(1:nd)'*x(1:nd) + c(nd+1:N)'*x(nd+1:N);
end

% MODIFICATION LOG:
%
% 060622 ango Wrote file
% 080603 med  Switched to conAssign, cleaned
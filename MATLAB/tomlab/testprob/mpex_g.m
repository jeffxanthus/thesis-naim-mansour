% function g = mpex_g(x, Prob)
%
% Compute gradient vector to problem Prob.P
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% g      Gradient g(x)

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Jun 22, 2006.  Last modified Jun 3, 2008.

function g = mpex_g(x,Prob)

P = Prob.P;
if P == 1
    % Gauvin and Savard
    g = [2*x(1); 2*(x(2)-10) ];
elseif P == 2
    % DF1
    g0 = x(1)-1-x(2);
    g = 2*[ g0 ; -g0 ];
elseif P == 3
    % BILIN
    g = [ 8 4 -4 40 4 0 0 0 ]';
elseif P == 4
    % Bard3m
    g = [ -2*x(1), -3, -4 , 2*x(4), 0, 0];
elseif P == 5
    % RalphMOD
    P = Prob.user.P;
    nd = 4;
    ns = 100;
    N = ns+nd;
    g = Prob.user.c;
    g(1:nd)   = g(1:nd) + ...
        P(1:nd,1:nd)*x(1:nd) + ...
        P(1:nd,nd+1:N)*x(nd+1:N);
    g(nd+1:N) = g(nd+1:N) + ...
        (x(1:nd)'*P(1:nd,nd+1:N))' + ...
        P(nd+1:N,nd+1:N)*x(nd+1:N);
end

% MODIFICATION LOG:
%
% 060622 ango Wrote file
% 080603 med  Switched to conAssign, cleaned
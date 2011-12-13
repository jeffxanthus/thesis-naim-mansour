% function dc = mpex_dc(x, Prob)
%
% Compute nonlinear constraints Jacobian to problem Prob.P
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% dc     Constraint Jacobian dc(x)

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Jun 22, 2006.  Last modified Jun 3, 2008.

function dc = mpex_dc(x,Prob)

P = Prob.P;
if P == 1
    % Gauvin and Savard
    dc = [ 4 8 1 ; -1 -1 0 ];
elseif P == 2
    % DF1
    dc = [ 2*x(1) , 0 ; 2*(x(1)-1) , 2*(x(2)-1) ; -2*x(1) , 1 ];
elseif P == 3
    % BILIN
    dc = [];
elseif P == 4
    % Bard3m
    dc = [ 2*x(1), 2, 0,0,0,0 ;  ...
        2*x(1)-2 ,   2*x(2), -2, 1, 0,0  ; ...
        0, 1, 3, -4, 0, 0 ; ...
        0, 0, 2,  0, 2, -3 ; ...
        0, 0, 0,  0, -1, 4 ; ...
        ];
elseif P == 5
    % RalphMOD
    dc = [];
end

% MODIFICATION LOG:
%
% 060622 ango Wrote file
% 080603 med  Switched to conAssign, cleaned
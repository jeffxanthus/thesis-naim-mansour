% function c = mpex_c(x, Prob)
%
% Compute nonlinear constraints to problem Prob.P
%
% x      Point x where f(x) is evaluated
% Prob   Problem structure
% c      Constraint vector c(x)
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Jun 22, 2006.  Last modified Jun 3, 2008.

function c = mpex_c(x,Prob)

P = Prob.P;
if P == 1
    % Gauvin and Savard
    %    Fy: 0 <= 4 * (x + 2*y - 30) + u   complements   y >= 0;
    %
    %    Fu: 0 <=  20 - x - y              complements   u >= 0;
    c = [ ...
        4 * (x(1) + 2*x(2) - 30) + x(3) ; ...
        20 - x(1) - x(2) ; ...
        ];
elseif P == 2
    % DF1
    c = [  x(1)^2  ; ...
        (x(1)-1)^2 + (x(2)-1)^2  ; ...
        x(2) - x(1)^2 + 1   ;  ...
        ];
elseif P == 3
    % BILIN
    c = [];
elseif P == 4
    % Bard3m
    c = [ ...
        x(1)^2 + 2*x(2) ; ...
        x(1)^2 - 2*x(1) +   x(2)^2 - 2*x(3) + x(4) + 3  ; ... %  _|_  x(5) >= 0;
        x(2)   + 3*x(3) - 4*x(4) - 4                  ; ... %  _|_  x(6) >= 0;
        (2*x(3)+2*x(5))-3*x(6) ; ... % _|_  x(3) >= 0;
        (-5-x(5))+4*x(6)       ; ... % _|_  x(4) >= 0;
        ];
elseif P == 5
    % RalphMOD
    c = [];
end

% MODIFICATION LOG:
%
% 060622 ango Wrote file
% 080603 med  Switched to conAssign, cleaned
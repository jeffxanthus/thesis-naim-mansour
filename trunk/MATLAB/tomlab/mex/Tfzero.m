% Tfzero, TOMLAB fzero routine
%
% Find a zero for f(x) in the interval [x_L, x_U]. Tfzero searches for a
% zero of a function f(x) between the given scalar values x_L and x_U until
% the width of the interval (xLow, xUpp) has collapsed to within a
% tolerance specified by the stopping criterion,
%     abs(xLow-xUpp) <= 2.*(RelErr*abs(xLow)+AbsErr).
% The method used is an efficient combination of bisection and the secant
% rule and is due to T. J. Dekker.
%
% Syntax of Tfzero:
%
% [xLow, xUpp, ExitFlag] = Tfzero(x_L, x_U, Prob, x_0, RelErr, AbsErr)
%
% INPUT:
% x_L         Lower limit on the zero x to f(x)
% x_U         Upper limit on the zero x to f(x)
% Prob        Structure, sent to Matlab routine ZeroFunc. The function name
%             should be set in Prob.FUNCS.f0. Only the function will be
%             used, not the gradient.
% x_0         An initial guess on the zero to f(x). If empty, x_0 is set as
%             the middle point in [x_L, x_U].
% RelErr      Relative error tolerance, default 1E-7
% AbsErr      Absolute error tolerance, default 1E-14
%
% OUTPUT:
% xLow        Lower limit on the zero x to f(x)
% xUpp        Upper limit on the zero x to f(x)
% ExitFlag    Status flag 1,2,3,4,5
%
%          1  xLow is within the requested tolerance of a zero. The
%             interval (xLow, xUpp) collapsed to the requested tolerance,
%             the function changes sign in (xLow, xUpp), and f(x) decreased
%             in magnitude as (xLow, xUpp) collapsed.
%
%          2  f(xLow) = 0. However, the interval (xLow, xUpp) may not have
%             collapsed to the requested tolerance.
%
%          3  xLow may be near a singular point of f(x). The interval
%             (xLow, xUpp) collapsed to the requested tolerance and the
%             function changes sign in (xLow, xUpp), but f(x) increased in
%             magnitude as (xLow, xUpp) collapsed, i.e. abs(f(xLow)) >
%             max(abs(f(xLow-IN)),abs(f(xUpp-IN))).
%
%          4  No change in sign of f(x) was found although the interval
%             (xLow, xUpp) collapsed to the requested tolerance. The user
%             must examine this case and decide whether xLow is near a
%             local minimum of f(x), or xLow is near a zero of even
%             multiplicity, or neither of these.
%
%          5  Too many (> 500) function evaluations used.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Sept 30, 2000.  Last modified Oct 10, 2007.
%# mex

function ...
    [xLow, xUpp, ExitFlag] = Tfzero(x_L, x_U, Prob, x_0, RelErr, AbsErr)
help Tfzero.m

% MODIFICATION LOG:
%
% 041201 med  Notation changed to x_L and x_U
% 041201 med  Added help that Prob.FUNCS.f0 should be set
% 060814 med  FUNCS used for callbacks instead
% 071010 med  Help updated, wich gradient
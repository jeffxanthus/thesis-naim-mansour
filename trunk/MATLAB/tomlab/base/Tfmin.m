% Tfmin, TOMLAB fminbnd routine
%
% Minimize function of one variable. Find miniumum x in [x_L, x_U] for
% function Func within tolerance xTol. Solves using Brents minimization
% algorithm.
%
% Syntax of Tfzero:
%
% [x, nFunc] = Tfmin(Func, x_L, x_U, xTol, Prob)
%
% INPUT:
% Func        Function of x to be minimized. Func must be defined as:
%             f = Func(x) if no 5th argument Prob is given
%             or
%             f = Func(x, Prob) if 5th argument Prob is given.
%
% x_L         Lower bound on x
% x_U         Upper bound on x
% xTol        Tolerance on accuracy of minimum
% Prob        Structure (or any Matlab variable) sent to Func. If many
%             parameters are to be sent to Func set them in Prob as a
%             structure. Example for parameters a and b:
%
%             Prob.user.a = a;
%             Prob.user.b = b;
%             [x, nFunc] = Tfmin('myFunc',0,1,1E-5,Prob);
%
%             In myFunc:
%
%             function f = myFunc(x, Prob)
%             a = Prob.user.a;
%             b = Prob.user.b;
%             f = "matlab expression dependent of x, a and b";
%
% OUTPUT:
% x           Solution
% nFunc       Number of calls to Func
%
% Reference: "Computer Methods for Mathematical Computations", Forsythe,
% Malcolm, and Moler, Prentice-Hall, 1976.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Nov 26, 2000.   Last modified Oct 10, 2007.

function [x, nFunc] = Tfmin(Func, x_L, x_U, xTol, Prob)

if nargin < 4
    [x, nFunc] = tomsol(8, x_L, x_U, Func, []);
else
    [x, nFunc] = tomsol(8, x_L, x_U, Func, xTol, Prob);
end

% MODIFICATION LOG:
%
% 030116 hkh  Change name to use Tfmin
% 041201 med  Notation changed to x_L and x_U
% 060714 med  Spelling updated
% 071010 med  Help updated
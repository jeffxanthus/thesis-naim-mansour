% geno_d2c.m
%
% function d2c = geno_d2c(x, Prob)
%
% Evaluates the second part of the Lagrangian for test functions from the
% GENO development manual.
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function d2c = geno_d2c(x, lam, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        d2c = [];
    case 2 % GENO The Colville #4 Function
        d2c = [];
    case 3 % GENO The 2-Dimensional Rastrigin Function
        d2c = [];
    case 4 % GENO The Economic Dispatch Problem
        d2c = [];
    case 5 % GENO A Pressure Vessel Design Problem
        d2c = zeros(4,4);
        d2c(3,3) = lam*(-pi*2*x(4)-4/3*pi*3*2*x(3));
        d2c(3,4) = lam*(-pi*2*x(3));
        d2c(4,3) = d2c(3,4);
        d2c = d2c/1000;
    case 6 % GENO The Alkylation Process
        d2c = [];
    case 7 % GENO Decentralised Economic Planning
        d2c = [];
    case 8 % GENO Heat Exchanger Optimisation
        d2c = [];
    case 9 % GENO The Harvest Problem
        d2c = [];
    case 10 % GENO A Non-linear Resource Allocation Problem
        d2c = [];
    case 11 % GENO Oligopolist Market Equilibrium Problem
        d2c = [];
    case 12 % GENO The Euclidean Compromise Solution I
        d2c = [];
    case 13 % GENO The Euclidean Compromise Solution II
        d2c = [];
    case 14 % GENO Efficient Portfolio Selection
        d2c = [];
    case 15 % GENO A Dynamic Non-cooperative Game 1
        d2c = [];
    case 16 % GENO A Dynamic Non-cooperative Game 2
        d2c = [];
    case 17 % GENO A Dynamic Non-cooperative Game 3
        d2c = [];
    case 18 % GENO A Dynamic Non-cooperative Game 4
        d2c = [];
    case 19 % GENO A Dynamic Non-cooperative Game 5
        d2c = [];
    case 20 % GENO A Dynamic Non-cooperative Game 6
        d2c = [];
    case 21 % GENO A Dynamic Non-cooperative Game 7
        d2c = [];
end

% MODIFICATION LOG
%
% 060802  med  Created, based on GENO manual
% 080603  med  Switched to *Assign, cleaned
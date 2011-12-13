% geno_dc.m
%
% function dc = geno_dc(x, Prob)
%
% Evaluates the constraint Jacobian for test functions from the GENO
% development manual.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function dc = geno_dc(x, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        dc = zeros(3,5);
        dc(1,1) = 0.0006262*x(4);
        dc(1,2) = 0.0056858*x(5);
        dc(1,3) = -0.0022053*x(5);
        dc(1,4) = 0.0006262*x(1);
        dc(1,5) = -0.0022053*x(3) + 0.0056858*x(2);
        dc(2,1) = 0.0029955*x(2);
        dc(2,2) = 0.0071317*x(5) + 0.0029955*x(1);
        dc(2,3) = 2*0.0021813*x(3);
        dc(2,5) = 0.0071317*x(2);
        dc(3,1) = 0.0012547*x(3);
        dc(3,3) = 0.0047026*x(5) + 0.0012547*x(1) + 0.0019085*x(4);
        dc(3,4) = 0.0019085*x(3);
        dc(3,5) = 0.0047026*x(3);
    case 2 % GENO The Colville #4 Function
        dc = [];
    case 3 % GENO The 2-Dimensional Rastrigin Function
        dc = [];
    case 4 % GENO The Economic Dispatch Problem
        dc = zeros(3,4);
        dc(1,1) = -1;
        dc(1,3) = -1000*cos(-x(3) - 0.25);
        dc(1,4) = -1000*cos(-x(4) - 0.25);
        dc(2,2) = -1;
        dc(2,3) = 1000*cos(x(3) - 0.25) + 1000*cos(x(3) -x(4) - 0.25);
        dc(2,4) = -1000*cos(x(3) -x(4) - 0.25);
        dc(3,3) = -1000*cos(x(4) -x(3) - 0.25);
        dc(3,4) = 1000*cos(x(4) - 0.25) + 1000*cos(x(4) -x(3) - 0.25);
    case 5 % GENO A Pressure Vessel Design Problem
        dc = zeros(1,4);
        dc(1,3) = -pi*2*x(3)*x(4)-4/3*pi*3*x(3)^2;
        dc(1,4) = -pi*x(3)^2;
        dc = dc/1000;
    case 6 % GENO The Alkylation Process
        dc = zeros(6,10);
        dc(1,1) = 1.12 + 0.13167*x(8) - 0.00667*(x(8)^2);
        dc(1,4) = -0.99;
        dc(1,8) = 0.13167*x(1) - 2*0.00667*x(1)*x(8);
        dc(2,6) = 0.325;
        dc(2,7) = -0.99;
        dc(2,8) = 1.098 - 2*0.038*x(8);
        dc(3,1) = - dc(1,1);
        dc(3,4) = ((1.0/0.99) - 0.99) - dc(1,4);
        dc(3,8) = -dc(1,8);
        dc(4,6) = - dc(2,6);
        dc(4,7) = ((1.0/0.99) - 0.99) - dc(2,7);
        dc(4,8) = -dc(2,8);
        dc(5,3) = 98000*(x(4)*x(9)+1000*x(3))^(-1) ...
            -98000*x(3)*(x(4)*x(9)+1000*x(3))^(-2)*1000;
        dc(5,4) = -98000*x(3)*(x(4)*x(9)+1000*x(3))^(-2)*x(9);
        dc(5,6) = -1;
        dc(5,9) = -98000*x(3)*(x(4)*x(9)+1000*x(3))^(-2)*x(4);
        dc(6,1) = -(x(2)+x(5))*(x(1))^(-2);
        dc(6,2) = 1/x(1);
        dc(6,5) = 1/x(1);
        dc(6,8) = -1;
    case 7 % GENO Decentralised Economic Planning
        dc = zeros(4,7);
        dc(1,1) = -4*x(1);
        dc(1,2) = -12*x(2)^3;
        dc(1,3) = -1;
        dc(1,4) = -8*x(4);
        dc(1,5) = -5;
        dc(2,1) = -7;
        dc(2,2) = -3;
        dc(2,3) = -20*x(3);
        dc(2,4) = -1;
        dc(2,5) = -1;
        dc(3,1) = -23;
        dc(3,2) = -2*x(2);
        dc(3,6) = -12*x(6);
        dc(3,7) = 8;
        dc(4,1) = -8*x(1)+3*x(2);
        dc(4,2) = -2*x(2)+3*x(1);
        dc(4,3) = -4*x(3);
        dc(4,6) = -5;
        dc(4,7) = 11;
    case 8 % GENO Heat Exchanger Optimisation
        dc = zeros(3,8);
        dc(1,1) = x(6)-100;
        dc(1,4) = -833.33252;
        dc(1,6) = x(1);
        dc(2,2) = x(7)-x(4);
        dc(2,4) = -x(2)+1250;
        dc(2,5) = -1250;
        dc(2,7) = x(2);
        dc(3,3) = x(8)-x(5);
        dc(3,5) = -x(3)+2500;
        dc(3,8) = x(3);
    case 9 % GENO The Harvest Problem
        dc = [];
    case 10 % GENO A Non-linear Resource Allocation Problem
        dc = [];
    case 11 % GENO Oligopolist Market Equilibrium Problem
        dc = [];
    case 12 % GENO The Euclidean Compromise Solution I
        dc = [];
    case 13 % GENO The Euclidean Compromise Solution II
        dc = [];
    case 14 % GENO Efficient Portfolio Selection
        dc = [];
    case 15 % GENO A Dynamic Non-cooperative Game 1
        dc = [];
    case 16 % GENO A Dynamic Non-cooperative Game 2
        dc = [];
    case 17 % GENO A Dynamic Non-cooperative Game 3
        dc = [];
    case 18 % GENO A Dynamic Non-cooperative Game 4
        dc = [];
    case 19 % GENO A Dynamic Non-cooperative Game 5
        dc = [];
    case 20 % GENO A Dynamic Non-cooperative Game 6
        T = 5;
        dc = zeros(T,11);
        for i=1:5
            dc(i,T+1+i) = dc(i,T+1+i) - 1;
            dc(i,T+i) = dc(i,T+i) + 1 + x(i) - 2*x(i) + 2*1.25*(x(T+i))*x(i) - 0.75*(x(T+i))^2*x(i);
            dc(i,i) = dc(i,i) + (x(T+i) - 2*x(T+i) + 1.25*(x(T+i))^2 - 0.25*(x(T+i))^3);
        end
    case 21 % GENO A Dynamic Non-cooperative Game 7
        dc = [];
end

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 080603  med  Switched to *Assign, cleaned
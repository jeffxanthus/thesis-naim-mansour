% geno_c.m
%
% function cx = geno_c(x, Prob)
%
% Evaluates the constraints for test functions from the GENO
% development manual.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function cx = geno_c(x, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        cx = zeros(3,1);
        cx(1,1) = 85.334407 + 0.0056858*x(2)*x(5) + 0.0006262*x(1)*x(4) - 0.0022053*x(3)*x(5);
        cx(2,1) = 80.512490 + 0.0071317*x(2)*x(5) + 0.0029955*x(1)*x(2) + 0.0021813*x(3)*x(3);
        cx(3,1) = 9.3009610 + 0.0047026*x(3)*x(5) + 0.0012547*x(1)*x(3) + 0.0019085*x(3)*x(4);
    case 2 % GENO The Colville #4 Function
        cx = [];
    case 3 % GENO The 2-Dimensional Rastrigin Function
        cx = [];
    case 4 % GENO The Economic Dispatch Problem
        cx = zeros(3,1);
        cx(1,1) = 1000*sin(-x(3) - 0.25) + 1000*sin(-x(4) - 0.25) + 894.8 - x(1);
        cx(2,1) = 1000*sin(x(3) - 0.25) + 1000*sin(x(3) -x(4) - 0.25) + 894.8 - x(2);
        cx(3,1) = 1000*sin(x(4) - 0.25) + 1000*sin(x(4) -x(3) - 0.25) + 1294.8;
    case 5 % GENO A Pressure Vessel Design Problem
        cx = -pi*x(3)^2*x(4)-4/3*pi*x(3)^3;
        cx = cx/1000;
    case 6 % GENO The Alkylation Process
        cx = zeros(6,1);
        cx(1,1) = 1.12*x(1) + 0.13167*x(1)*x(8) - 0.00667*x(1)*(x(8)^2) - 0.99*x(4);
        cx(2,1) = 57.425 + 1.098*x(8) - 0.038*(x(8)^2) + 0.325*x(6) - 0.99*x(7);
        cx(3,1) = ((1.0/0.99) - 0.99)*x(4) - cx(1,1);
        cx(4,1) = ((1.0/0.99) - 0.99)*x(7) - cx(2,1);
        cx(5,1) = -x(6)+98000*x(3)/(x(4)*x(9)+1000*x(3));
        cx(6,1) = -x(8) + (x(2)+x(5))/x(1);
    case 7 % GENO Decentralised Economic Planning
        cx = zeros(4,1);
        cx(1,1) = 127 - 2*x(1)^2 - 3*(x(2)^4) - 4*(x(4)^2) - 5*x(5) - x(3);
        cx(2,1) = 282 - 7*x(1) - 3*x(2) - 10*x(3)^2 - x(4) - x(5);
        cx(3,1) = 196 - 23*x(1) - x(2)^2 - 6*x(6)^2 + 8*x(7);
        cx(4,1) = -4*x(1)^2 - x(2)^2 + 3*x(1)*x(2) - 2*x(3)^2 + 11*x(7) - 5*x(6);
    case 8 % GENO Heat Exchanger Optimisation
        cx = zeros(3,1);
        cx(1,1) = x(1)*x(6) - 833.33252*x(4) - 100*x(1);
        cx(2,1) = x(2)*x(7) - 1250*x(5) - x(2)*x(4) + 1250*x(4);
        cx(3,1) = x(3)*x(8) - x(3)*x(5) + 2500*x(5);
    case 9 % GENO The Harvest Problem
        cx = [];
    case 10 % GENO A Non-linear Resource Allocation Problem
        cx = [];
    case 11 % GENO Oligopolist Market Equilibrium Problem
        ci = [10;8;6;4;2];
        Ki = [5;5;5;5;5];
        alphai = [1/1.2; 1/1.1; 1.00; 1/0.9; 1/0.8];
        beta = 1/1.1;
        Q = sum(x);
        pQ = (5000/Q).^(beta);
        fq = ci.*x+(x.^(1+alphai)./(1+alphai)./Ki.^alphai);
        c1 = x.*pQ-fq;
        cx = [c1(5)-c1(4); c1(4)-c1(3); c1(3)-c1(2); c1(2)-c1(1)];
    case 12 % GENO The Euclidean Compromise Solution I
        cx = [];
    case 13 % GENO The Euclidean Compromise Solution II
        cx = [];
    case 14 % GENO Efficient Portfolio Selection
        cx = [];
    case 15 % GENO A Dynamic Non-cooperative Game 1
        cx = [];
    case 16 % GENO A Dynamic Non-cooperative Game 2
        cx = [];
    case 17 % GENO A Dynamic Non-cooperative Game 3
        cx = [];
    case 18 % GENO A Dynamic Non-cooperative Game 4
        cx = [];
    case 19 % GENO A Dynamic Non-cooperative Game 5
        cx = [];
    case 20 % GENO A Dynamic Non-cooperative Game 6
        T = 5;
        cx = zeros(T,1);
        for i=1:5
            cx(i,1) = -x(T+1+i) + x(T+i) + (x(T+i) - 2*x(T+i) + 1.25*(x(T+i))^2 - 0.25*(x(T+i))^3)*x(i);
        end
    case 21 % GENO A Dynamic Non-cooperative Game 7
        y  = x(4:7); % Integer variables
        cx = [x(1)^2+x(2)^2+x(3)^2 + y(3)^2; ...
            x(2)^2+y(2)^2; x(3)^2+y(3)^2; x(3)^2+y(2)^2];
end

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 060807  med  Corrected heat exchanger problem
% 080603  med  Switched to *Assign, cleaned
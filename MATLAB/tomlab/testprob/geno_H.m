% geno_H.m
%
% function H = geno_H(x, Prob)
%
% Evaluates the Hessian for test functions from the GENO development
% manual.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function H = geno_H(x, Prob)

x=x(:);
P=Prob.P;

switch P
    case 1 % GENO Ex 1
        H = [];
    case 2 % GENO The Colville #4 Function
        H = [];
    case 3 % GENO The 2-Dimensional Rastrigin Function
        H = [];
    case 4 % GENO The Economic Dispatch Problem
        H = [];
    case 5 % GENO A Pressure Vessel Design Problem
        H = zeros(4,4);
        H(1,1) = 3.1661*0.0625^2*x(4)*2 + 19.84*0.0625^2*x(3)*2;
        H(1,3) = 0.6224*0.0625*x(4) + 19.84*0.0625^2*2*x(1);
        H(3,1) = H(1,3);
        H(1,4) = 0.6224*0.0625*x(3) + 3.1661*0.0625^2*2*x(1);
        H(4,1) = H(1,4);
        H(2,3) = 2*0.0625*1.7781*x(3);
        H(3,2) = H(2,3);
        H(3,4) = 0.6224*x(1)*0.0625;
        H(4,3) = H(3,4);
        H(3,3) = 1.7781*0.0625*x(2)*2;
    case 6 % GENO The Alkylation Process
        H = [];
    case 7 % GENO Decentralised Economic Planning
        H = [];
    case 8 % GENO Heat Exchanger Optimisation
        H = [];
    case 9 % GENO The Harvest Problem
        H = [];
    case 10 % GENO A Non-linear Resource Allocation Problem
        H = zeros(11,11);
        for j=1:5
            for i=1:5
                if i~=j
                    H(j,i) = 1;
                    for k=1:5
                        if k~=i & k~=j
                            H(j,i) = H(j,i)*(1+k*x(k));
                        else
                            H(j,i) = H(j,i)*k;
                        end
                    end
                else
                    H(j,i) = 0;
                end
            end
        end
        H = -H;
    case 11 % GENO Oligopolist Market Equilibrium Problem
        H = [];
    case 12 % GENO The Euclidean Compromise Solution I
        H = [];
    case 13 % GENO The Euclidean Compromise Solution II
        H = [];
    case 14 % GENO Efficient Portfolio Selection
        H = [];
    case 15 % GENO A Dynamic Non-cooperative Game 1
        H = [];
    case 16 % GENO A Dynamic Non-cooperative Game 2
        H = [];
    case 17 % GENO A Dynamic Non-cooperative Game 3
        H = [];
    case 18 % GENO A Dynamic Non-cooperative Game 4
        H = [];
    case 19 % GENO A Dynamic Non-cooperative Game 5
        H = [];
    case 20 % GENO A Dynamic Non-cooperative Game 6
        H = [];
    case 21 % GENO A Dynamic Non-cooperative Game 7
        H = [];
end

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 080603  med  Switched to *Assign, cleaned
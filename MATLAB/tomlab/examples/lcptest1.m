% lcptest1 
% Transportation problem in TOMLAB 
%
% 0 <= x_k _|_ Mx_k+q >= 0  (_|_ means complements)
%
% SUPPLY                  DEMAND
% Seattle:   325          New York:  325
% San Deigo: 575          Chicago:   300
%                         Topeka:    275
%
% TRANSPORT COSTS:
%            New York   Chicago   Topeka
% Seattle    0.225      0.153     0.162
% San Deigo  0.225      0.162     0.126
%
% The variabels x_k are identified as:
% x(1)  = Amount shipped from Seattle   -> New York
% x(2)  = Amount shipped from Seattle   -> Chicago
% x(3)  = Amount shipped from Seattle   -> Topeka
% x(4)  = Amount shipped from San Diego -> New York
% x(5)  = Amount shipped from San Diego -> Chicago
% x(6)  = Amount shipped from San Diego -> Topeka
% 
% x(7)  = Price at market New York
% x(8)  = Price at market Chicago
% x(9)  = Price at market Topeka
% x(10) = Price at plant Seattle
% x(11) = Price at plant San Diego

% clear all

M = [ 0  0  0  0  0  0 -1  0  0  1  0; % Seattle   -> New York
      0  0  0  0  0  0  0 -1  0  1  0; % Seattle   -> Chicago
      0  0  0  0  0  0  0  0 -1  1  0; % Seattle   -> Topeka
      0  0  0  0  0  0 -1  0  0  0  1; % San Deigo -> New York
      0  0  0  0  0  0  0 -1  0  0  1; % San Diego -> Chicago
      0  0  0  0  0  0  0  0 -1  0  1; % San Diego -> Topeka
      1  0  0  1  0  0  0  0  0  0  0; 
      0  1  0  0  1  0  0  0  0  0  0; 
      0  0  1  0  0  1  0  0  0  0  0; 
     -1 -1 -1  0  0  0  0  0  0  0  0; 
      0  0  0 -1 -1 -1  0  0  0  0  0];

% The q vector shows the costs and demand/supply.

q = [ 0.225; 0.153; 0.162; 0.225; 0.162; 0.126; -325.000; -300.000; -275.000; 325.000; 575.000 ];

% Meaning of first complimentary condition (first row in M with q):
% 1. If New York price - Seattle price + cost == 0, then x_k(1) may be
% anywhere between bounds. Anything quantity be shipped.
% 2. If New York price - Seattle price + cost > 0, then x_k(1) must be
% zero (or low lomit). Ship minimum to make the most money.
% 3. If New York price - Seattle price + cost < 0, then x_k(1) must be at
% upper bounds. Ship maximum to make the most money.

% Meaning of seventh complimentary condition (seventh row in M with q):
% 1. If amount shipped from Seattle to New York + amount shipped from San
% Diego to New York - Demand in New York == 0, then x_k(7) may be anywhere
% between bounds.
% 2. If amount shipped from Seattle to New York + amount shipped from San
% Diego to New York - Demand in New York > 0, then x_k(7) must be
% zero (or low lomit). I.e. free in a rational market.
% 3. If amount shipped from Seattle to New York + amount shipped from San
% Diego to New York - Demand in New York < 0, then x_k(7) must be at
% upper bounds. Increase price to make more money.

% Similar for other parts.

% Assign the problem.

c = q;

x_L = zeros(11,1);
x_U = inf*ones(11,1);

A   = M;
b_L = -q;
b_U = inf*ones(length(q),1);

MPEC = [1 0 1 0 0 0; ...
        2 0 2 0 0 0; 
        3 0 3 0 0 0; 
        4 0 4 0 0 0; 
        5 0 5 0 0 0; 
        6 0 6 0 0 0; 
        7 0 7 0 0 0; 
        8 0 8 0 0 0; 
        9 0 9 0 0 0; 
        10 0 10 0 0 0;         
        11 0 11 0 0 0];

Name = 'lcpTest1';

x_0 = [];

Prob = lcpAssign(c, x_L, x_U, x_0, A, b_L, b_U, MPEC, Name);

% Solve the problem.
Prob.KNITRO.options.ALG = 3;
Prob.PriLevOpt = 2;
Result = tomRun('knitro', Prob, 2);

% The first 11 original variables
x = Result.x_k(1:11)

% The linear rows in M*x+q
Mxq = Prob.A(:,1:11) * x + q

% Complementarity? This should be (close to) zero: 

x'*Mxq


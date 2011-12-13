% geno_prob.m
%
% Defines constrained global optimization problems from GENO manual.
%
% function [probList, Prob] = geno_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function [probList, Prob] = geno_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'GENO Ex 1'...
    ,'GENO The Colville #4 Function'...
    ,'GENO The 2-Dimensional Rastrigin Function'...
    ,'GENO The Economic Dispatch Problem'...
    ,'GENO A Pressure Vessel Design Problem'...
    ,'GENO The Alkylation Process'...
    ,'GENO Decentralised Economic Planning'...
    ,'GENO Heat Exchanger Optimisation'...
    ,'GENO The Harvest Problem'...
    ,'GENO A Non-linear Resource Allocation Problem'...
    ,'GENO Oligopolist Market Equilibrium Problem'...
    ,'GENO The Euclidean Compromise Solution I'...
    ,'GENO The Euclidean Compromise Solution II'...
    ,'GENO Efficient Portfolio Selection'...
    ,'GENO A Dynamic Non-cooperative Game 1'...
    ,'GENO A Dynamic Non-cooperative Game 2'...
    ,'GENO A Dynamic Non-cooperative Game 3'...
    ,'GENO A Dynamic Non-cooperative Game 4'...
    ,'GENO A Dynamic Non-cooperative Game 5'...
    ,'GENO A Dynamic Non-cooperative Game 6'...
    ,'GENO A Dynamic Non-cooperative Game 7'); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_0 = [];
IntVars = [];
isLS = 0;

if P == 1
    Name = 'GENO Ex 1';
    x_L = [78 33 27 27 27]';
    x_U = [102 45 45 45 45]';
    % Variable 2 shoule be restricted to chunks of 0.34
    b_L = [];
    b_U = [];
    A = [];
    c_L = [0;90;20];
    c_U = [92;110;25];
    x_opt = [0 0 0 0 0];
    f_opt = 0;
    x_min = [-1 -1];
    x_max = [ 1  1];
elseif P == 2
    Name = 'GENO The Colville #4 Function';
    x_L = -10*ones(4,1);
    x_U =  10*ones(4,1);
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = [];
    x_opt = [1 1 1 1];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 3
    Name = 'GENO The 2-Dimensional Rastrigin Function';
    x_L = [-1 -1]';
    x_U = [ 1  1]';
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = [];
    x_opt = [0 0];
    f_opt = -2.0;
    x_min = x_L;
    x_max = x_U;
elseif P == 4
    Name = 'GENO The Economic Dispatch Problem';
    x_L = [0 1e-6 -0.55 -0.55]';
    x_U = [1200 1200 0.55 0.55]';
    b_L = [-0.55; -0.55];
    b_U = [inf inf];
    A = [0 0 -1 1; 0 0 1 -1];
    c_L = [0 0 0]';
    c_U = [0 0 0]';
    x_opt = [589.265856 1124.772836 0.184406 -0.365594];
    f_opt = 4221.956525;
    x_min = x_L;
    x_max = x_U;
    x_0 = 0.1*ones(4,1);
elseif P == 5
    Name = 'GENO A Pressure Vessel Design Problem';
    x_L = [1 1 10 10]';
    x12_U = round(99/0.0625);
    x_U = [x12_U x12_U 200 200]';
    % x(4) - 240 <= 0 redundant constraint
    b_L = [-inf; -inf];
    b_U = [0; 0];
    A = [-0.0625 0 0.0193 0; ...
        0 -0.0625 0.00954 0];
    c_L = -inf;
    c_U = -1296;
    x_opt = [13 7 42.098446 176.636596];
    f_opt = 6059.714336;
    IntVars = logical([1;1;0;0]);
    x_min = x_L;
    x_max = x_U;
elseif P == 6
    Name = 'GENO The Alkylation Process';
    x_L = [0.00001 0.00001 0.00001 0.00001 0.00001 85 90 3 1.2 145]';
    x_U = [2000 16000 120 5000 2000 93 95 12 4 162]';
    b_L = [-35.82 133 35.82 -133 0];
    b_U = [inf inf inf inf 0];
    A = [0 0 0 0 0 0 0 0 -0.9 -0.222; ...
        0 0 0 0 0 0 3 0 0 -0.99; ...
        0 0 0 0 0 0 0 0 1/0.9-0.9+0.9 0.222; ...
        0 0 0 0 0 0 -3 0 0 1/0.99-0.99+0.99; ...
        -1 0 0 1.22 -1 0 0 0 0 0];
    c_L = [0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;0;0];
    x_opt = [1695.747560 15766.851262 54.144966 3027.278310 ...
        1997.531978 90.084472 94.984529 10.475842 1.571577 153.485638];
    f_opt = -1763.745454;
    x_min = x_L;
    x_max = x_U;
    x_0 = 0.5*(x_U+x_L);
elseif P == 7
    Name = 'GENO Decentralised Economic Planning';
    x_L = -10*ones(7,1);
    x_U = 10*ones(7,1);
    b_L = []; b_U = []; A = [];
    c_L = [0;0;0;0];
    c_U = [inf;inf;inf;inf];
    x_opt = [2.330500 1.951372 -0.477541 4.365727 -0.624487 1.038131 1.594227];
    f_opt = 680.630057;
    x_min = x_L;
    x_max = x_U;
elseif P == 8
    Name = 'GENO Heat Exchanger Optimisation';
    x_L = [100 1000 1000 10 10 10 10 10]';
    x_U = [10000 10000 10000 1000 1000 1000 1000 1000]';
    b_L = [-1 -1 -1]';
    b_U = [inf inf inf]';
    A = [0 0 0 -0.0025 0 -0.0025 0 0; ...
        0 0 0 0.0025 -0.0025 0 -0.0025 0; ...
        0 0 0 0 0.01 0 0 -0.01];
    c_L = [-83333.333;0;1250000];
    c_U = [inf;inf;inf];
    x_opt = [579.306665 1359.970680 5109.970675 182.017720 295.601180 ...
        217.982302 286.416525 395.601173];
    f_opt = 7049.248021;
    x_min = x_L;
    x_max = x_U;
elseif P == 9
    Name = 'GENO The Harvest Problem';
    x_L = [zeros(4,1);100;-inf;-inf;-inf;100]';
    x_U = [inf*ones(4,1);100;inf;inf;inf;100]';
    b_L = zeros(4,1);
    b_U = zeros(4,1);
    A = zeros(4,9);
    for i=1:4
        A(i,[i i+4 i+5]) = [-1 1.1 -1];
    end
    c_L = [];
    c_U = [];
    x_opt = [7.513148 9.090909 11 13.31 100 102.486852 103.644628 103.009091 100];
    f_opt = -12.721038;
    x_min = x_L;
    x_max = x_U;
    x_0 = 0.1*ones(9,1);
elseif P == 10
    Name = 'GENO A Non-linear Resource Allocation Problem';
    x_L = [0 0 0 0 0 8 -inf -inf -inf -inf -inf]';
    x_U = [5 5 5 5 5 8  inf  inf  inf  inf  inf]';
    b_L = [zeros(5,1);8];
    b_U = [zeros(5,1);8];
    A = zeros(6,11);
    for i=1:5
        A(i,[i i+5 i+6]) = [-1 1 -1];
    end
    A(end,1:5) = ones(1,5);
    c_L = [];
    c_U = [];
    x_opt = [1 1 2 2 2 8 7 6 4 2 0];
    f_opt = -4158.000000;
    IntVars = ones(11,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 11
    Name = 'GENO Oligopolist Market Equilibrium Problem';
    x_L = ones(5,1);
    x_U = 100*ones(5,1);
    b_L = []; b_U = []; A = [];
    c_L = zeros(4,1);
    c_U = inf*ones(4,1);
    x_opt = [36.932511 41.818141 43.706578 42.659240 39.178952];
    f_opt = [];
    x_min = x_L;
    x_max = x_U;
    x_0 = [36.1;63.5;40.5;40.6;21.7]';
    y   = zeros(5,1);
    t   = zeros(5,1);
    isLS = 1;
elseif P == 12
    Name = 'GENO The Euclidean Compromise Solution I';
    x_L = -6;
    x_U = 6;
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = [];
    x_opt = 1;
    f_opt = [1 1]';
    x_min = x_L;
    x_max = x_U;
    t = zeros(2,1);
    y = zeros(2,1);
    isLS = 1;
elseif P == 13
    Name = 'GENO The Euclidean Compromise Solution II';
    x_L = [-5 -5]';
    x_U = [ 5  5]';
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = [];
    x_opt = [2.5 2.5];
    f_opt = 0.5*(2.5^2*2);
    x_min = x_L;
    x_max = x_U;
    t = zeros(2,1);
    y = zeros(2,1);
    isLS = 1;
elseif P == 14
    Name = 'GENO Efficient Portfolio Selection';
    x_L = zeros(6+7,1);
    x_U = ones(6+7,1);
    b_L = [zeros(6,1);1];
    b_U = [zeros(6,1);1];
    A   = zeros(7,13);
    for i=1:6
        A(i,[i i+6 i+7]) = [-1 1 -1];
    end
    A(end,1:6) = ones(1,6);
    c_L = [];
    c_U = [];
    x_opt = [0 0 0.3 0 0.5 0.2 1 1 1 0.7 0.7 0.2 0];
    f_opt = [];
    x_min = x_L;
    x_max = x_U;
elseif P == 15
    Name = 'GENO A Dynamic Non-cooperative Game 1';
    x_L = [-100*ones(10,1);0;-inf*ones(5,1);0;-inf*ones(5,1)]';
    x_U = [100*ones(10,1);0;inf*ones(5,1);0;inf*ones(5,1)]';
    b_L = zeros(10,1);
    b_U = zeros(10,1);
    A = zeros(10,22);
    for i=1:5
        A(i,[11+i 16+i]) = [-1 1];
        A(5+i, [17+i 16+i 10+i 5+i]) = [-1 2 -1 1/25];
    end
    c_L = [];
    c_U = [];
    x_opt = [0 0 0 0 0 0.8 0.6 0.4 0.2 0 0 0 0.032 0.088 0.16 0.24 ...
        0 0.032 0.088 0.16 0.24 0.32];
    f_opt = -0.120000;
    x_min = x_L;
    x_max = x_U;
elseif P == 16
    Name = 'GENO A Dynamic Non-cooperative Game 2';
    x_L = [-1*ones(10,1);1;-2*ones(5,1);2;-2*ones(5,1)];
    x_U = [ones(10,1);1;2*ones(5,1);2;2*ones(5,1)];
    b_L = zeros(10,1);
    b_U = zeros(10,1);
    A = zeros(10,22);
    for i=1:5
        A(i,[11+i 10+i 16+i i]) = [-1 1 1 1];
        A(i+5,[17+i 16+i 10+i 5+i]) = [-1 1 -1 1];
    end
    c_L = [];
    c_U = [];
    x_opt = [-1 -1 0 1 0 -1 1 1 1 0 1 2 1 0 0 0 2 0 -1 -1 0 0];
    f_opt = 11.5;
    x_min = x_L;
    x_max = x_U;
elseif P == 17
    Name = 'GENO A Dynamic Non-cooperative Game 3';
    x_L = [zeros(5,1);1.0;zeros(5,1)];
    x_U = [inf*ones(5,1);1;inf*ones(5,1)];
    b_L = [zeros(5,1);-inf*ones(5,1)];
    b_U = zeros(10,1);
    A = zeros(10,11);
    for i=1:5
        A(i,[i+6 i+5 i]) = [-1 1.02 -1.02];
        A(i+5, [i i+5]) = [1 -1];
    end
    c_L = [];
    c_U = [];
    x_opt = [0.225661 0.211887 0.198954 0.186810 0.175407 ...
        1 0.789825 0.589497 0.398354 0.215775 0.041175];
    f_opt = -2.105094;
    x_min = x_L;
    x_max = x_U;
    x_0 = ones(11,1)*0.1;
elseif P == 18
    Name = 'GENO A Dynamic Non-cooperative Game 4';
    x_L = [-inf*ones(5,1);100;-inf*ones(5,1)];
    x_U = [inf*ones(5,1);100;inf*ones(5,1)];
    b_L = zeros(5,1);
    b_U = zeros(5,1);
    A = zeros(5,11);
    for i=1:5
        A(i,[i+6 i+5 i]) = [-1 0.01 1];
    end
    c_L = [];
    c_U = [];
    x_opt = [-0.500013 -0.0025 -0.000012 0 0 100 0.499987 0.002499 0 0 0];
    f_opt = 10000.500012;
    x_min = x_L;
    x_max = x_U;
elseif P == 19
    Name = 'GENO A Dynamic Non-cooperative Game 5';
    x_L = [ones(5,1);8;-inf*ones(5,1)];
    x_U = [5*ones(5,1);8;inf*ones(5,1)];
    b_L = [zeros(5,1);8];
    b_U = [zeros(5,1);8];
    A = zeros(6,11);
    for i=1:5
        A(i,[i+6 i+5 i]) = [-1 1 -1];
    end
    A(end,1:5) = ones(1,5);
    c_L = [];
    c_U = [];
    x_opt = [1.056631 1.556632 1.723301 1.806636 1.856800 8 6.943369 5.386737 3.663436 1.856800 0];
    f_opt = -4415.712521;
    x_min = x_L;
    x_max = x_U;
elseif P == 20
    Name = 'GENO A Dynamic Non-cooperative Game 6';
    x_L = [-1*ones(5,1);2;zeros(5,1)];
    x_U = [1*ones(5,1);2;3*ones(5,1)];
    A = [];
    b_L = [];
    b_U = [];
    c_L = zeros(5,1);
    c_U = zeros(5,1);
    x_opt = [0.963664 0.022570 0.001489 -0.571663 -1.005097 2 2.963664 2.997695 2.999929 2.142424 1];
    f_opt = 0.629013;
    x_min = x_L;
    x_max = x_U;
elseif P == 21
    Name = 'GENO A Dynamic Non-cooperative Game 7';
    x_L = [zeros(5,1);100;-inf*ones(4,1);100];
    x_U = [inf*ones(5,1);100;inf*ones(4,1);100];
    b_L = zeros(5,1);
    b_U = zeros(5,1);
    A = zeros(5,11);
    for i=1:5
        A(i,[i+6 i+5 i]) = [-1 1.1 -1];
    end
    c_L = [];
    c_U = [];
    x_opt = [6.830134 8.264463 10 12.1 14.641 100 103.169866 105.222389 105.744628 104.219091 100];
    f_opt = -15.955390;
    x_min = x_L;
    x_max = x_U;
    x_0 = ones(11,1);
else
    error('geno_prob: Illegal problem number');
end

% Set x_0 to zeros (dummy for GUI)
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the routines to compute the function f and the constraints c.
if P == 14
    g = [];
else
    g = 'geno_g';
end

c = 'geno_c';
dc = 'geno_dc';

if isempty(c_L) & isempty(c_U)
    c = [];
    dc = [];
end

if (P == 5 | P == 10)
    H = 'geno_H';
    d2c = 'geno_d2c';
else
    H = [];
    d2c = [];
end

if isLS
    Prob = clsAssign('geno_f', 'geno_g', [], x_L, x_U, Name, x_0, ...
        y, t, [], [], [], [], ...
        A, b_L, b_U, c, dc, [], c_L, c_U, ...
        x_min, x_max, f_opt, x_opt, IntVars);
else
    Prob = minlpAssign('geno_f', g, H, [], x_L, x_U, Name, x_0, ...
        IntVars, [], [], [], ...
        A, b_L, b_U, c, dc, d2c, [], c_L, c_U,...
        x_min, x_max, f_opt, x_opt);
end
Prob.P = P;

% MODIFICATION LOG
%
% 060801  med  Created, based on GENO manual
% 070221  hkh  Use logical to set IntVars with 0/1
% 080603  med  Switched to *Assign, cleaned
% 080916  nhq  Switched all x_opt into row-vectors

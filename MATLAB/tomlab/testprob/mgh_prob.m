% mgh_prob
%
% Defines test problems from More, Garbow, Hillstrom nonlinear least squares.
%
% The following can be change in the file:
%
% uP(1)    Scale the initial x_0, x_0=uP(1)*x_0. uP(1) in [-100000,100000]
% uP(2)    Problem dimension n, set for problem 20 - 35, i.e.
%          '20. Watson [More G H #20]' to ,'35. Chebyquad'
% uP(3)    Number of residuals m, set for problem 32 - 35, i.e.
%          '32. Linear - Full Rank' to '35. Chebyquad'
%
% If hard coded variable is changed to: mgh_bounds = 1;
% then simple bound constraints are added following D. GAY in
% his report "A TRUST-REGION APPROACH TO LINEARLY
% CONSTRAINED OPTIMIZATION" for the following problems
% 1, 2, 6, 7, 8, 10, 12, 13, 15, 16, 17, 19, 20, 27, 32, 33, 34, 35
% DEFAULT: mgh_bounds = 0;
%
% function [probList, Prob] = mgh_prob(P, ask, Prob);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function [probList, Prob] = mgh_prob(P, varargin)

mgh_bounds = 0; % Flag if adding bounds

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Rosenbrock'...
    ,'Freudenstein and Roth'...
    ,'Powell Badly Scaled'...
    ,'Brown Badly Scaled'...
    ,'Beale'...
    ,'Jennrich and Sampson'...
    ,'Helical Valley'...
    ,'Bard'...
    ,'Gaussian'...
    ,'Meyer'...
    ,'Gulf Research and Development'...
    ,'Box 3-Dimensional'...
    ,'Powell Singular'...
    ,'Wood'...
    ,'Kowalik and Osborne'...
    ,'Brown and Dennis'...
    ,'Osborne1'...
    ,'Biggs'...
    ,'Osborne2'...
    ,'Watson'...
    ,'Extended Rosenbrock'...
    ,'Extended Powell Singular'...
    ,'Penalty I'...
    ,'Penalty II'...
    ,'Variably Dimensioned'...
    ,'Trigonometric'...
    ,'Brown Almost Linear'...
    ,'Discrete Boundary Value'...
    ,'Discrete Integral Equation'...
    ,'Broyden Tridiagonal'...
    ,'Broyden Banded'...
    ,'Linear - Full Rank'...
    ,'Linear - Rank 1'...
    ,'Linear - Rank 1 w 0 cols & rows'...
    ,'Chebyquad'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

JacPattern = []; t = []; weightType = []; weightY = [];
SepAlg = []; f_Low = []; A = []; b_L = []; b_U = [];
uP = [];

if P> 20 & P < 36
    [Name, t, y, x_0, x_opt, f_opt, x_L, x_U, x_min, x_max, ...
        JacPattern, uP] = mgh2(P, ProbDef, mgh_bounds);
elseif P==1 % Simple bounds added
    Name='Rosenbrock [More G H #1]';
    % x_0 Scale Factor
    x0S = 1;
    y = [0 0]';
    x_0 = [-1.2 1]';
    x_0 = x0S*x_0;
    if mgh_bounds
        x_opt = [0.5 0.25];
        f_opt = 0.125;
    else
        x_opt = [1 1];
        f_opt = 0;
    end
    x_L = [-50 0]';
    x_U = [0.5 100]';
    x_min = [-1.5 -2]';
    x_max = [ 1.5  2]';
elseif P==2 % Simple bounds added
    Name='Freudenstein and Roth [More G H #2]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [11.37 -0.9];
        f_opt = 24.492601;
    else
        x_opt = [5 4; ...
            11.412779 -0.896805];
        f_opt = [0; 24.49212684];
    end
    y = [0 0]';
    x_0 = [0.5 -2]';
    x_0 = x0S*x_0;
    x_L = [0 -30]';
    x_U = [20 -0.9]';
    x_min = [0 -3]';
    x_max = [6  5]';
elseif P==3
    Name='Powell Badly Scaled [More G H #3]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds used
    x_opt = [1.098e-5 9.106];
    f_opt = 0;
    y = [0 0]';
    x_0 = [0 1]';
    x_0 = x0S*x_0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [0 -3]';
    x_max = [1e-4  10]';
elseif P==4
    Name='Brown Badly Scaled [More G H #4]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds used
    x_0 = [1 1]';
    x_0 = x0S*x_0;
    x_opt = [1e6 2e-6];
    y = [0 0 0]';
    f_opt = 0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [0 0]';
    x_max = [1e7 1]';
elseif P==5
    Name='Beale [More G H #5]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds used
    y = [-1.5 -2.25 -2.625]';
    x_0 = [1 1]';
    x_0 = x0S*x_0;
    x_opt = [3 0.5];
    f_opt = 0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [0 0]';
    x_max = [4 2]';
elseif P==6 % Simple bounds added
    Name='Jennrich and Sampson [More G H #6]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [0.26000 0.255751];
        f_opt = 62.191556458188948;
    else
        x_opt = [0.257825 0.257825];
        f_opt = 62.181091178;
    end
    y = zeros(10,1);
    x_0 = [0.3 0.4]';
    x_0 = x0S*x_0;
    x_L = [0.26 0]';
    x_U = [10 20 ]';
    x_min = [0 0 ]';
    x_max = [0.5 0.5]';
elseif P==7 % Simple bounds added
    Name='Helical Valley [More G H #7]';
    % x_0 Scale Factor
    x0S = 1;
    y = zeros(3,1);
    if mgh_bounds
        f_opt = 0.495211060480699850;
        x_opt = [0.8 0.561937 0.964934];
    else
        f_opt = 0;
        x_opt = [1 0 0];
    end
    x_0 = [-1 0 0]';
    x_0 = x0S*x_0;
    x_L = [  -100 -1 -1 ]';
    x_U = [   0.8  1  1 ]';
    x_min = [-2 -1 -1 ]';
    x_max = [2 1 1]';
elseif P==8 % Simple bounds added
    Name='Bard [More G H #8]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [0.1;1.519451;1.981874];
        f_opt = 0.00479114236;
    else
        x_opt = [0.082411;1.133036;2.343695];
        f_opt = 0.5*8.21487e-3;
    end
    y = [0.14 0.18 0.22 0.25 0.29 0.32 0.35 0.39 0.37 0.58 0.73 0.96 1.34 2.10 4.39]';
    x_0 = [1 1 1]';
    x_0 = x0S*x_0;
    x_L = [0.1 0 0]';
    x_U = [50 100 50]';
    x_min = [0 0 0]';
    x_max = [2 2 2]';
elseif P==9
    Name='Gaussian [More G H #9]';
    % x_0 Scale Factor
    x0S = 1;
    y = [0.0009 0.0044 0.0175 0.0540 0.1295 0.2420 0.3521 0.3989 0.3521 ...
        0.2420 0.1295 0.0540 0.0175 0.0044 0.0009]';
    m=length(y);
    t = (8-(1:m)')./2;
    x_0 = [0.4 1 0]';
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0.5*1.12793e-8;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [0 0 0]';
    x_max = [2 2 2]';
elseif P==10 % Simple bounds added
    Name='Meyer [More G H #10]';
    % x_0 Scale Factor
    if mgh_bounds
        x_opt = [0.006;6125.418985;343.337723];
        f_opt = 63.7258856;
    else
        x_opt = [0.005610;6181.34629;345.223633];
        f_opt = 0.5*8.79458e1;
    end
    y = [34780 28610 23650 19630 16370 13720 11540  9744 ...
        8261  7030  6005  5147  4427  3820  3307  2872]';
    m=length(y);
    x_0 = [0.02 4000 250]';
    t = 45+5*(1:m)';
    x_L = [0.006 0 0]';
    x_U = [2 3E5 3E4]';
    x_min = [0 3000 200]';
    x_max = [0.1 5000 300]';
elseif P==11
    Name='Gulf Research and Development [More G H #11]';
    % x_0 Scale Factor
    x0S = 1;
    m=10;
    y = zeros(m,1);
    t = 0.01*(1:m)';
    % Sligthly changed starting point
    x_0 = [5 2.5 0.15]'; % clsSolve converges to other numerical stat.point
    x_0 = x0S*x_0;
    x_opt = [50 25 1.5];
    f_opt = 0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [0 0 0]';
    x_max = [60 30 2]';
elseif P==12 % Simple bounds added
    Name='Box 3-Dimensional [More G H #12]';
    % x_0 Scale Factor
    x0S = 1;
    m=10;
    if mgh_bounds
        x_opt = [1.037944;9.5;0.972228];
        f_opt = 0.00005717409;
    else
        x_opt = [1 10 1;10 1 -1;1 1 0];
        f_opt = zeros(3,1);
    end
    t = 0.1*(1:m)';
    y = zeros(m,1);
    x_0 = [0 10 20]';
    x_0 = x0S*x_0;
    x_L = [0 5 0]';
    x_U = [2 9.5 20]';
    x_min = [0 0 0]';
    x_max = [50 50 50]';
elseif P==13 % Simple bounds added
    Name='Powell Singular [More G H #13]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [0.1;-0.009982;0.043073;0.043784];
        f_opt = 0.000093909815;
    else
        x_opt = zeros(1,4);
        f_opt = 0;
    end
    y = zeros(4,1);
    x_0 = [3 -1 0 1]';
    x_0 = x0S*x_0;
    x_L = [0.1 -20 -1 -1]';
    x_U = [100 20 1 50]';
    x_min = [-1 -1 -1 -1]';
    x_max = [5 5 5 5]';
elseif P==14
    Name='Wood [More G H #14]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds
    y = zeros(6,1);
    x_0 = [-3 -1 -3 -1]';
    x_0 = x0S*x_0;
    x_opt = [1 1 1 1];
    f_opt = 0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = [-4 -2 -4 -2]';
    x_max = [2 2 2 2]';
elseif P==15 % Simple bounds added
    Name='Kowalik and Osborne [More G H #15]';
    % x_0 Scale Factor
    x0S = 1;
    % n = 4, m = 11
    if mgh_bounds
        x_opt = [0.193033 0.197840 0.13 0.138346];
        f_opt = 0.000153915361732099;
    else
        x_opt = [0.192807 0.191282 0.123057 0.136062;Inf -14.07 -Inf -Inf];
        f_opt = 0.5*[3.07505e-4;1.02734E-3];
    end
    y = [0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246]';
    x_0 = [0.25 0.39 0.415 0.39]';
    x_0 = x0S*x_0;
    x_L = [0 -1 0.13 0.12]';
    x_U = [10 12 13 12]';
    x_min = [0 0 0 0]';
    x_max = [1 1 1 1]';
elseif P==16 % Simple bounds added
    Name='Brown and Dennis [More G H #16]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [-10;12.677818;-0.485528;0.2];
        f_opt = 44430.239883750226;
    else
        x_opt = [-11.594440;13.203630;-0.403439;0.236779];
        f_opt = 0.5*8.58222e4;
    end
    m=20;
    t = 0.2*(1:m)';
    y = zeros(m,1);
    x_0 = [25 5 -5 -1]';
    x_0 = x0S*x_0;
    x_L = [-10  0 -100 -20  ]';
    x_U = [100 15    0   0.2]';
    x_min = -ones(length(x_0),1);
    x_max =  ones(length(x_0),1);
elseif P==17 % Simple bounds added
    Name='Osborne1 [More G H #17]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [0.375092 1.9 -1.4286 0.012793 0.022272];
        f_opt = 0.000027351416716868;
    else
        x_opt = [0.3754 1.9358 -1.4647 0.01287 0.02212];
        f_opt = 0.5*5.46489e-5;
    end
    y = [0.844 0.908 0.932 0.936 0.925 0.908 0.881 0.850 0.818 0.784 0.751 0.718 ...
        0.685 0.658 0.628 0.603 0.580 0.558 0.538 0.522 0.506 0.490 0.478 0.467 ...
        0.457 0.448 0.438 0.431 0.424 0.420 0.414 0.411 0.406]';
    m=length(y);
    t = 10*((1:m)-1)';
    x_0 = [0.5 1.5 -1 0.01 0.02]';
    x_0 = x0S*x_0;
    x_L = [ 0 0   -50  0   0 ]';
    x_U = [50 1.9 -0.1 10 10 ]';
    x_min = -ones(length(x_0),1);
    x_max =  ones(length(x_0),1);
elseif P==18
    Name='Biggs [More G H #18]';
    % x_0 Scale Factor
    x0S = 1;
    m=13;
    t = 0.1*(1:m)';
    y = exp(-t)+3.*exp(-4.*t)-5.*exp(-10.*t);
    x_0 = [1 2 1 1 1 1]';
    x_0 = x0S*x_0;
    x_opt = [1 10 1 5 4 3];
    f_opt = 0;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    x_min = -ones(length(x_0),1);
    x_max =  ones(length(x_0),1);
elseif P==19 % Simple bounds added
    Name='Osborne2 [More G H #19]';
    % x_0 Scale Factor
    x0S = 1;
    if mgh_bounds
        x_opt = [1 0.5 0.684698 0.6 0.901119 0.340623 0.712022...
            123.324574 2.169916 5.108346 0];
        f_opt = 0.255191627055202760;
    else
        x_opt = [1.3100 0.4315 0.6336 0.5993 0.7539 0.9056 1.3651 ...
            4.8248 2.3988 4.5689 5.6754];
        f_opt = 0.5*4.01683e-2;
    end
    y = [1.366 1.191 1.112 1.013 0.991 0.885 0.831 0.847 0.786 0.725 0.746 0.679 ...
        0.608 0.655 0.615 0.606 0.602 0.626 0.651 0.724 0.649 0.649 0.694 0.644 ...
        0.624 0.661 0.612 0.558 0.533 0.495 0.500 0.423 0.395 0.375 0.372 0.391 ...
        0.396 0.405 0.428 0.429 0.523 0.562 0.607 0.653 0.672 0.708 0.633 0.668 ...
        0.645 0.632 0.591 0.559 0.597 0.625 0.739 0.710 0.729 0.720 0.636 0.581 ...
        0.428 0.292 0.162 0.098 0.054]';
    m=length(y);
    t = 0.1*((1:m)-1)';
    x_0 = [1.3 0.65 0.65 0.7 0.6 3 5 7 2 4.5 5.5]';
    x_0 = x0S*x_0;
    n = length(x_0);
    x_L = zeros(n,1);
    x_U = 100*ones(n,1);
    x_L(1)=1; x_L(2)=0.5; x_L(4)=0.6;
    x_U(1)=150; x_U(6)=500; x_U(7)=500;
    % NOTE THAT x_U(11)==0, so the starting point is infeasible
    x_U(8)=500; x_U(10)=10; x_U(11)=0;
    x_min = -ones(length(x_0),1);
    x_max =  ones(length(x_0),1);
elseif P==20 % Simple bounds added
    Name='Watson [More G H #20]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable 2 <= n <= 31, default n=6
    m=31;
    n = 20; % Give problem dimension (2-31) (often 6,9,12 or 20)
    uP(1)=n;
    uP(2)=n;
    t = (1:m)'./29;
    y = zeros(m,1);
    x_0 = zeros(n,1);
    x_0 = x0S*x_0;
    if mgh_bounds
        x_L = -10*ones(n,1);
        x_L(1) = -0.015;
        x_U = 100*ones(n,1);
        x_U(1) = 10;
        x_opt = [];
        if n==6
            % Solution from clsSolve. Slightly better than MINOS and NPSOL.
            f_opt = 0.001144803219124709;
            x_U(6) = 0.99;
        elseif n==9
            % Solution from MINOS.
            f_opt = 0.018700698896592250;
            x_L = [-1E-5  0    0    0  0 -3  0 -3  0 ]';
            x_U = [ 1E-5  0.9  0.1  1  1  0  4  0  2 ]';
        elseif n==12
            % Solution from MINOS.
            x_L = [-1  0 -1 -1   -1  0 -3  0 -10  0 -5  0 ]';
            x_U = [ 0  0  0  0.3  0  1  0 10   0 10  0  1 ]';
            x_opt = [-0.00601;0;0;0.3;0;1;...
                0;4.901949;-10;3.703968;0;0.728223];
            f_opt = 5.8971220901104857;
        else
            f_opt = [];
        end
    else
        x_opt = [];
        x_L=[]; x_U=[];
        if n==6
            f_opt = 0.5*2.28767e-3;

        elseif n==9
            f_opt = 0.5*1.39976e-6;
        elseif n==12
            x_opt = [0;1.000002;-0.000564;0.347821;-0.156732;1.052815;...
                -3.247271;7.288435;-10.271848;9.074114;-4.541375;1.012012];
            f_opt = 0.5*4.72238e-10;
        else
            f_opt = [];
        end
    end
    x_min = -ones(n,1);
    x_max =  ones(n,1);
else
    error('mgh_prob: Illegal problem number');
end

if ~mgh_bounds
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
end

Prob = clsAssign('mgh_r','mgh_J', JacPattern, x_L, x_U, Name, x_0, ...
    y, t, weightType, weightY, SepAlg, f_Low, ...
    A, b_L, b_U, [], [], [], [], [], ...
    x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.uP = uP;

function [Name, t, y, x_0, x_opt, f_opt, x_L, x_U, x_min, x_max, ...
    JacPattern, uP] = mgh2(P, Prob, mgh_bounds)

t=[];
uP=[];
JacPattern=[];

if P==21
    Name='Extended Rosenbrock [More G H #21]';
    uP=checkuP(Name,Prob);
    % x_0 Scale Factor
    x0S = 1;
    % No bounds
    % n variable but even, 2 <= n , default n=30
    n = 30;
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    x_0 = zeros(n,1);
    x_0(1:2:n-1) = -1.2;
    x_0(2:2:n)   =  1;
    x_0 = x0S*x_0;
    x_opt = ones(n,1);
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==22
    Name='Extended Powell Singular [More G H #22]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds
    % n variable but multiple of 4, 4 <= n , default n=40
    n = 40; % Give problem dimension (multiple of 4);
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    x_0 = zeros(n,1);
    x_0(1:4:n) =  3;
    x_0(2:4:n) = -1;
    x_0(3:4:n) =  0;
    x_0(4:4:n) =  1;
    x_0 = x0S*x_0;
    x_opt = zeros(n,1);
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==23
    Name='Penalty I [More G H #23]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n+1;
    y = zeros(m,1);
    x_0 = (1:n)';
    x_0 = x0S*x_0;
    x_opt = [];
    if n==4
        f_opt = 0.5*2.24997e-5;
    elseif n==10
        f_opt = 0.5*7.08765e-5;
    else
        f_opt = [];
    end
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==24
    Name='Penalty II [More G H #24]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds
    % n variable, 1 <= n, default n=30
    n = 30; % Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=2*n;
    y = zeros(m,1);
    x_0 = 0.5*ones(n,1);
    x_0 = x0S*x_0;
    x_opt = [];
    if n==4
        f_opt = 0.5*9.37629e-6;
    elseif n==10
        f_opt = 0.5*2.93660e-4;
    else
        f_opt = [];
    end
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==25
    Name='Variably Dimensioned [More G H #25]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(2)=n;
    m=n+2;
    y = zeros(m,1);
    x_0 = 1-(1:n)'./n;
    x_0 = x0S*x_0;
    x_opt = ones(n,1);
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==26
    Name='Trigonometric [More G H #26]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    x_0   = 1/n*ones(n,1);
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==27 % Simple bounds added
    Name='Brown Almost Linear [More G H #27]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    if mgh_bounds
        if n==10
            f_opt = 0.003046589288298371;
            x_opt = [0.969674;1;0.9;0.969674*ones(6,1);1.348716];
        elseif n==12
            f_opt = 0.003408552998621086;
            x_opt = [0.977037;1;0.9;0.977037*ones(8,1);1.336157];
        else
            f_opt = [];
            x_opt = [];
        end
    else
        x_opt = ones(n,1);
        f_opt = 0;
    end
    y = zeros(m,1);
    x_0   = 0.5*ones(n,1);
    x_0 = x0S*x_0;
    x_L = zeros(n,1);
    x_U = 100*ones(n,1);
    if n >= 3
        x_L(2)=1.0; x_U(3)=0.9;
        if n > 10
            for i = 11:n
                x_L(i) = x_L(i-10);
                x_U(i) = x_U(i-10);
            end
        end
    end
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==28
    Name='Discrete Boundary Value [More G H #28]';
    % x_0 Scale Factor
    x0S = 1;
    % No bounds
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    t = (1:n)'./(n+1);
    x_0 = t.*( t - 1);
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==29
    Name='Discrete Integral Equation [More G H #29]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    t = (1:n)'./(n+1);
    x_0 = t.*( t - 1);
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==30
    Name='Broyden Tridiagonal [More G H #30]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    x_0 = -ones(n,1);
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==31
    Name='Broyden Banded [More G H #31]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m=n;
    y = zeros(m,1);
    x_0 = -ones(n,1);
    x_0 = x0S*x_0;
    x_opt = [];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==32 % Simple bounds added
    Name='Linear - Full Rank [More G H #32]';
    % n variable, 1 <= n, default n=30
    % m variable, n <= m, default m=50
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m = 50; %Give m (m >= n)
    uP(3)=m;
    if mgh_bounds
        x_opt = [-0.5;-ones(n-1,1)];
        f_opt = 0.125+0.5*(m-n);
    else
        x_opt = -ones(n,1);
        f_opt = 0.5*(m-n);
    end
    y = zeros(m,1);
    x_0 = ones(n,1);
    x_L = [-0.5; -2*ones(n-1,1)];
    x_U = 100*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==33 % Simple bounds added
    Name='Linear - Rank 1 [More G H #33]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    % m variable, n <= m, default m=50
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m = 50; %Give m (m >= n)
    uP(3)=m;
    if mgh_bounds
        % constrained minimum seems the same as the unconstrained
        x_opt = [];
        f_opt = 0.5*m*(m-1)/(2*(2*m+1));
    else
        x_opt = [];
        f_opt = 0.5*m*(m-1)/(2*(2*m+1));
    end
    y = zeros(m,1);
    x_0 = ones(n,1);
    x_0 = x0S*x_0;
    x_L = [-0.9; -2*ones(n-1,1)];
    x_U = 100*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==34 % Simple bounds added
    Name='Linear - Rank 1 w 0 cols & rows [More G H #34]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=30
    % m variable, n <= m, default m=50
    n = 30; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m = 50; %Give m (m >= n)
    uP(3)=m;
    if mgh_bounds
        x_opt = [];
        f_opt = [];
        x_L = zeros(n,1);
        x_L(2)=-0.4;
        x_L(3)=0.3;
        x_U = 100*ones(n,1);
    else
        x_L=[]; x_U=[];
        x_opt = [];
        f_opt = 0.5*(m^2+3*m-6)/(2*(2*m-3));
    end
    y = zeros(m,1);
    x_0 = ones(n,1);
    x_0 = x0S*x_0;
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==35 % Simple bounds added
    Name='Chebyquad [More G H #35]';
    % x_0 Scale Factor
    x0S = 1;
    % n variable, 1 <= n, default n=8
    % m variable, n <= m, default m=8
    n = 8; %Give problem dimension
    uP(1)=n;
    uP(2)=n;
    m = 8; %Give m (m >= n)
    uP(3)=m;
    y = zeros(m,1);
    x_0 = (1:n)'.*(n+1);
    x_0 = x0S*x_0;
    if n == 1
        x_L = 0.5; x_U = 100;
    else
        x_L = zeros(n,1); x_U = ones(n,1);
        if n >= 3
            x_L(3)=0.1; x_U(1)=0.04; x_U(2)=0.2; x_U(3)=0.3;
            if (n >= 4) & (n ~= 8)
                x_U(1)=1.0; x_U(3)=0.23; x_U(4)=0.4;
                if n >= 10
                    x_L(2)=0.1; x_L(3)=0.2;  x_L(6)=0.5; x_L(7)=0.5; x_L(8)=0.5;
                    x_L(9)=0.5; x_L(10)=0.5; x_U(3)=0.3; x_U(5)=0.4;
                end
            end
        end
    end
    if mgh_bounds
        if n==m & n==8
            f_opt = 0.001796121192987531;
            x_opt = [0.04;0.188771;0.265004;0.731951;0.955861;0.497424;...
                0.497424;0.804466];
        else
            f_opt = [];
            x_opt = [];
        end
    else
        if (m==n) & (n<=10)
            if ((1<=n) & (n<=7)) | n==9
                f_opt = 0;
            elseif (n==8)
                f_opt = 0.5*3.51687e-3;
            elseif (n==10)
                f_opt = 0.5*6.50395e-3;
            end
        else
            f_opt = [];
        end
        x_opt = [];
    end
    x_min = -ones(n,1);
    x_max =  ones(n,1);
end

% MODIFICATION LOG:
%
% 990623  hkh  Use clsVarDef and clsProbSet instead of ProbVarDef and ProbSet
% 040414  hkh  First test if to use mgh2
% 041117  med  xxx_prob removed and code added
% 050303  hkh  uP(1) scaling of x_0. Add comments about uP settings
% 080603  med  Switched to clsAssign, cleaned
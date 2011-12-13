% glc_prob.m
%
% Defines constrained global optimization problems.
%
% function [probList, Prob] = glc_prob(P);
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
% Written June 1, 1999.   Last modified Sept 29, 2008.

function [probList, Prob] = glc_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
     'Gomez 2'...
    ,'Gomez 3'...
    ,'Hock-Schittkowski 59'...
    ,'Hock-Schittkowski 65'...
    ,'Hock-Schittkowski 104'...
    ,'Hock-Schittkowski 105'...
    ,'Schittkowski 234'...
    ,'Schittkowski 236'...
    ,'Schittkowski 237'...
    ,'Schittkowski 239'...
    ,'Schittkowski 330'...
    ,'Schittkowski 332'...
    ,'Schittkowski 343'...
    ,'Floudas-Pardalos 3.2 TP 1'...
    ,'Floudas-Pardalos 3.3 TP 2'...
    ,'Floudas-Pardalos 3.5 TP 4'...
    ,'Floudas-Pardalos 4.10 TP 9'...
    ,'Zimmerman'...
    ,'Bump 2'...
    ,'Bump 10'...
    ,'Bump 20'...
    ,'HGO 468:1 + constraint'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES
    
if isempty(P)
    return;
end

IntVars = []; f_Low = [];
x_0  = []; % Currently not used
user = [];
if P == 1
    Name = 'Gomez 2';
    x_L = [-1 -1]';
    x_U = [ 1  1]';
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = 0;
    x_opt = [0 0];
    f_opt = 0;
    x_min = [-1 -1];
    x_max = [ 1  1];
elseif P == 2
    Name = 'Gomez 3';
    x_L = [-1 -1]';
    x_U = [ 1  1]';
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = 0;
    x_opt = [0.109136094502324 -0.623423878335696];
    f_opt = -0.971103892157052;
    x_min = [-1 -1];
    x_max = [ 1  1];
elseif P == 3
    Name = 'Hock-Schittkowski 59';
    u = [75.196     3.8112    0.0020567  1.0345E-5   6.8306  0.030234   1.28134E-3 ...
        2.266E-7  0.25645   0.0034604  1.3514E-5  28.106   5.2375E-6  6.3E-8     ...
        7E-10     3.405E-4  1.6638E-6  2.8673     3.5256E-5];
    user.u = u;
    x_L = [0 0]';
    x_U = [75 65]';
    b_L = []; b_U = []; A = [];
    c_L = [0 0 0];
    c_U = [];
%     x_opt = [13.55010424 51.66018129];
%     f_opt = -7.804226324;
    x_opt = [13.551245575708313  51.655785872749753];
    f_opt = -7.802785527345220;
    x_min = x_L;
    x_max = x_U;
    x_0 = [90 10]'; % If running local solver
elseif P == 4
    Name = 'Hock-Schittkowski 65';
    x_L = [-4.5 -4.5 -5]';
    x_U = [ 4.5  4.5  5]';
    b_L = []; b_U = []; A = [];
    c_L = 0;
    c_U = [];
    x_opt = [3.650461821,3.650461821,4.6204170507];
    f_opt = 0.9535288567;
    x_min = x_L;
    x_max = x_U;
    x_0 = [-5 5 0]';
elseif P == 5
    Name = 'Hock-Schittkowski 104';
    % Orginal bounds - DIRECT method will fail to find feasible point in
    % reasonable time frame
    x_L = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';
    x_U = [7 3 1 1 7 7 2 1]';
    b_L = []; b_U = []; A = [];
    c_L = [0 0 0 0 1 -4.2];
    c_U = [];
    x_opt = [6.465114 2.232709 0.6673975 0.5957564 ...
        5.932676 5.527235 1.013322  0.4006682];
    f_opt = 3.9511634396;
    x_min = x_L;
    x_max = x_U;
    x_0 = [6 3 0.4 0.2 6 6 1 0.5]';
elseif P == 6
    Name = 'Hock-Schittkowski 105';
    y  = zeros(235,1);
    y(1)=95; y(2)=105; y(3:6)=110; y(7:10)=115; y(11:25)=120;
    y(26:40)=125;  y(41:55)=130; y(56:68)=135; y(69:89)=140;
    y(90:101)=145; y(102:118)=150; y(119:122)=155; y(123:142)=160;
    y(143:150)=165; y(151:167)=170; y(168:175)=175; y(176:181)=180;
    y(182:187)=185; y(188:194)=190; y(195:198)=195; y(199:201)=200;
    y(202:204)=205; y(205:212)=210; y(213)=215; y(214:219)=220;
    y(220:224)=230; y(225)=235; y(226:232)=240; y(233)=245;
    y(234:235)=250;
    user.y=y;
    x_L = [0.001 0.001 100 130 170 5 5 5]';
    x_U = [0.499 0.499 180 210 240 25 25 25]';
    b_L = -inf; b_U = 1;
    A = [1 1 0 0 0 0 0 0];
    c_L = []; c_U = [];
    x_opt = [0.4128928 0.4033526 131.2613 164.3135 217.4222 ...
        12.28018 15.7717 20.74682];
    f_opt = 1136;
    x_min = x_L;
    x_max = x_U;
    x_0 = [0.1 0.2 100 125 175 11.2 13.2 15.8]';
elseif P == 7
    Name = 'Schittkowski 234';
    x_L = [0.2 0.2]';
    x_U = [ 2   2 ]';
    b_L = []; b_U = []; A = [];
    c_L = 0;
    c_U = [];
    x_opt = [0.2 0.2];
    f_opt = -0.8;
    x_min = x_L;
    x_max = x_U;
    x_0 = [0 0]';
elseif P == 8
    Name = 'Schittkowski 236';
    B = [75.1963666677  -3.8112755343   0.1269366345  -2.0567665E-3   1.0345E-5 ...
        -6.8306567613   3.02344793E-2 -1.2813448E-3   3.52559E-5    -2.266E-7 ...
        0.2564581253  -3.460403E-3    1.35139E-5   -28.1064434908  -5.2375E-6 ...
        -6.3E-9         7E-10          3.405462E-4   -1.6638E-6     -2.8673112392]';
    user.B = B;
    x_L = [0 0]';
    x_U = [75 65]';
    b_L = []; b_U = []; A = [];
    c_L = [0;0];
    c_U = [];
    x_opt = [75 65];
    f_opt = -58.9034;
    x_min = x_L;
    x_max = x_U;
    x_0 = [90 10]';
elseif P == 9
    Name = 'Schittkowski 237';
    B = [75.1963666677  -3.8112755343   0.1269366345  -2.0567665E-3   1.0345E-5 ...
        -6.8306567613   3.02344793E-2 -1.2813448E-3   3.52559E-5    -2.266E-7 ...
        0.2564581253  -3.460403E-3    1.35139E-5   -28.1064434908  -5.2375E-6 ...
        -6.3E-9         7E-10          3.405462E-4   -1.6638E-6     -2.8673112392]';
    user.B = B;
    x_L = [54 0]';
    x_U = [75 65]';
    b_L = []; b_U = []; A = [];
    c_L = [0;0;0];
    c_U = [];
    x_opt = [75 65];
    f_opt = -58.9034;
    x_min = x_L;
    x_max = x_U;
    x_0 = [95 10]';
elseif P == 10
    Name = 'Schittkowski 239';
    B = [75.1963666677  -3.8112755343   0.1269366345  -2.0567665E-3   1.0345E-5 ...
        -6.8306567613   3.02344793E-2 -1.2813448E-3   3.52559E-5    -2.266E-7 ...
        0.2564581253  -3.460403E-3    1.35139E-5   -28.1064434908  -5.2375E-6 ...
        -6.3E-9         7E-10          3.405462E-4   -1.6638E-6     -2.8673112392]';
    user.B = B;
    x_L = [0 0]';
    x_U = [75 65]';
    b_L = []; b_U = []; A = [];
    c_L = 0;
    c_U = [];
    x_opt = [75 65];
    f_opt = -58.9034;
    x_min = x_L;
    x_max = x_U;
    x_0 = [95 10]';
elseif P == 11
    Name = 'Schittkowski 330';
    % x_L = [0 0]';  CHANGE BOUNDS AWAY FROM 0 TO AVOID OBVIOUS DIVISION WITH 0
    x_L = [1E-10 1E-10]';
    x_U = [5 5]';
    b_L = []; b_U = []; A = [];
    c_L = 0;
    c_U = [];
    x_opt = [1.287 0.5305];
    f_opt = 1.62058;
    x_min = x_L;
    x_max = x_U;
    x_0 = [2.5 2.5]';
elseif P == 12
    Name = 'Schittkowski 332';
    tmp1 = (1:100)';
    tmp2 = pi*(1/3+(tmp1-1)/180);
    user.tmp2 = tmp2; % Used in glc_f and glc_c
    x_L = [0 0]';
    x_U = [1.5 1.5]';
    b_L = []; b_U = []; A = [];
    c_L = -30;
    c_U = 30;
    x_opt = [0.9114 0.02928];
    f_opt = 29.92437939227878;
    x_min = x_L;
    x_max = x_U;
    x_0 = [0.75 0.75]';
elseif P == 13
    Name = 'Schittkowski 343';
    x_L = [0 0 0]';
    x_U = [36 5 125]';
    b_L = []; b_U = []; A = [];
    c_L = [0;0];
    c_U = [];
    %NHQ  Infinitely many optimal solutions, parameterized by x1:
    %     0 < x1 <= 36 ,  x2 = 675/x1^2 ,  x3 = sqrt(4.19*1E6/x1^2)
    %     f(x1,x2,x3) =  -5.6847825
%     x_opt = [36 0.520833333333333 56.859693029052004 ;
%              35 0.551020408163265 58.484255687024920 ];

    f_opt = -5.6847825;
    for k = 1:36
       x_opt(k,:) = [k , 675/k^2 , sqrt(4.19*1E6/k^2)];
    end
    x_min = x_L;
    x_max = x_U;
    x_0 = [22.3 0.5 125]';
elseif P == 14
    Name = 'Floudas-Pardalos 3.2 TP 1';
    x_L = 10*[10 100 200 1    1   1   1   1]';
    x_U = 1000*[1 2 6 0.5 0.5 0.5 0.5 0.5]';
    A = [0 0 0  0.0025   0    0.0025    0    0
        0 0 0 -0.0025 0.0025    0   0.0025  0
        0 0 0     0   -0.01     0      0   0.01];
    b_L = -inf*ones(3,1);
    b_U = [1;1;1];
    c_L = [];
    c_U = [0;0;-1250000];
    x_opt = [579.31 1359.97 5109.97 182.02 295.6 217.98 286.42 395.60];
    f_opt = 7049.25;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 15
    Name = 'Floudas-Pardalos 3.3 TP 2';
    x_L = [78 33 27 27 27]';
    x_U = [102 45 45 45 45]';
    b_L = []; b_U = []; A = [];
    c_L = [-85.334407;9.48751;10.699039];
    c_U = [6.665593;29.48751;15.699039];
    x_opt = [78 33 29.9953 45 36.7758];
    f_opt = -30665.5387;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 16
    Name = 'Floudas-Pardalos 3.5 TP 4';
    x_L = [0 0 0]';
    x_U = [2 2 3]'; % Upper bounds on x(2) added by hkh, from lin.eq.2
    A = [1 1 1
        0 3 1];
    b_L = -inf*ones(2,1);
    b_U = [4 6]';
    bvtmp = [3 1 2]';
    rtmp  = [1.5 -0.5 -5]';
    c_L = 0.25*bvtmp'*bvtmp-rtmp'*rtmp;
    c_U = [];
    % Note! [2 0 0] is also a global minimum with same f(x)=-4, not in book
    x_opt = [0.5 0 3; 2 0 0];
    f_opt = [-4;-4];
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 17
    Name = 'Floudas-Pardalos 4.10 TP 9';
    x_L = [0 0]';
    x_U = [3 4]';
    A = []; b_L=[]; b_U=[];
    c_L = [-2 -36]';
    c_U = [];
    x_opt = [2.3295 3.1783];
    f_opt = -5.5079;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 18
    Name = 'Zimmerman';
    x_L = [0 0]';
    x_U = [100 100]';
    A   = [];
    b_L = [];
    b_U = [];
    c_L = [];
    c_U = [16 14]';
    x_opt = [7 2];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 19 | P == 20 | P == 21
    if P==19
        Name = 'Bump 2';
        N   = 2;
    elseif P == 20
        Name = 'Bump 10';
        N   = 10;
    elseif P == 21
        Name = 'Bump 20';
        N   = 20;
    end
    x_L = 1E-100*ones(N,1);
    x_U = 10*ones(N,1);
    A   = ones(1,N);
    b_L = -inf;
    b_U = 7.5*N;
    c_L = 0.75;
    c_U = [];
    x_opt = [];
    f_opt = -100;
    if N == 2
        x_opt = [1.600484343511321 0.468602347203170];
        f_opt = -0.364984535866739;
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 22
    Name  = 'HGO 468:1 + constraint';
    x_L   = [0 0]';
    x_U   = [1 1]';
    A     = [];
    b_L   = [];
    b_U   = [];
    x_0   = [];
    c_L   = [];
    c_U   = 0;
    x_opt = [0.732213924029196  0.646336093752646];
    f_opt = -1.825389961181424;
%     x_opt = [0.799784924384833 0.639646734068316];
%     f_opt = -2.011755710223728;
    f_Low = -5;
    x_min = [0 0];
    x_max = [1 1];
elseif P == 23
    Name = 'Schittkowski 237 altered';
    B = [75.1963666677  -3.8112755343   0.1269366345  -2.0567665E-3   1.0345E-5 ...
        -6.8306567613   3.02344793E-2 -1.2813448E-3   3.52559E-5    -2.266E-7 ...
        0.2564581253  -3.460403E-3    1.35139E-5   -28.1064434908  -5.2375E-6 ...
        -6.3E-9         7E-10          3.405462E-4   -1.6638E-6     -2.8673112392]';
    user.B = B;
    x_L = [0 0]';
    x_U = [100 100]';
    b_L = []; b_U = []; A = [];
    c_L = [0;0];
    c_U = [];
    x_opt = [78.879405239716149  64.073057968805543];
    f_opt = -60.038142572669969;
    x_min = x_L;
    x_max = x_U;
    x_0 = [95 10]';
else
    error('glc_prob: Illegal problem number');
end

% Set x_0 to zeros
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the Prob
c  = 'glc_c';
if isempty(c_L) & isempty(c_U)
    c  = [];
end

Prob = glcAssign('glc_f', x_L, x_U, Name, A, b_L, b_U, ...
    c, c_L, c_U, x_0, IntVars, [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.user = user;

% MODIFICATION LOG
%
% 990413  mbk  Moved problem 1 to glb_*, added Gomez 2
% 990416  mbk  Small changes in comments.
% 990427  mbk  Problem 'Zimmerman' added.
% 990623  hkh  Use glbVarDef and glbProbSet instead of ProbVarDef and ProbSet
% 020110  hkh  Change lower bounds to avoid obvious division with 0 for
%              P == 11, "Schittkowski 330"
% 020111  hkh  Change upper bounds for P==5, to get feasible in time
% 020323  hkh  Adding Floudas-Pardalos chapter 12 problems
% 031126  hkh  Initial x_0 set as x_L for integer variables
% 040306  hkh  Correction of IntVars to [0 1] in problem 24
% 041117  med  xxx_prob removed and code added
% 050405  hkh  More decimals in some optimal f(x)
% 080423  nhq  More decimals in x_opt an f_opt for glc_prob 2.
% 080603  med  Switched to glcAssign, cleaned
% 081917  nhq  Found new global optimum for P == 13, 'Schittkowski 343'.
% 081919  nhq  Added x_opt and f_opt for P == 17, 'Bump 2'.
% 080925  nhq  Improved x_opt and f_opt for problem 16.
% 080925  nhq  Removed all IP problems to newly created glcIP_prob.
%              Changed problem numbers as well to avoid empty spaces.
% 080929  hkh  Fixed incorrect definition of bump problem 19-21
% 081001  nhq  Changed f_opt for glc_prob 3, old value not correct.
% 090512  nhq  Changed f_opt and x_opt glc_prob 13, old values not correct.
%              Infinitely many solutions, parametrized by x1.

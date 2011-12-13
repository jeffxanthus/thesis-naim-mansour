% Defines global optimization problems with simple bounds.
%
% glb_prob: Defines global optimization problems with simple bounds.
%
% function [probList, Prob] = glb_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc. $Release: 7.4.0$
% Written June 1, 1999.   Last modified Dec 5, 2009.

function [probList, Prob] = glb_prob(P, varargin)

if nargin < 1
    P=[];
end

% HGO 468:1 means problem 1 page 468 in "Handbook of Global Optimization"

probList=str2mat(...
    'Shekel 5'...
    ,'Shekel 7'...
    ,'Shekel 10'...
    ,'Hartman 3'...
    ,'Hartman 6'...
    ,'Branin RCOS'...
    ,'Goldstein and Price'...
    ,'Six-Hump Camel'...
    ,'Two-Dimensional Shubert'...
    ,'Sphere model 2'...
    ,'Sphere model 5'...
    ,'Sphere model 10'...
    ,'Griewanks function 2'...
    ,'Griewanks function 5'...
    ,'Griewanks function 10'...
    ,'Shekels foxholes 2'...
    ,'Shekels foxholes 5'...
    ,'Shekels foxholes 10'...
    ,'Michalewiczs function 2'...
    ,'Michalewiczs function 5'...
    ,'Michalewiczs function 10'...
    ,'Langermans function 2'...
    ,'Langermans function 5'...
    ,'Langermans function 10'...
    ,'HGO 468:1'...
    ,'HGO 468:2'...
    ,'Hock-Schittkowski 5'...
    ,'RB BANANA'...
    ,'Myers smooth' ...
    ,'Myers smoothly fluctuating' ...
    ,'Myers strongly fluctuating' ...
    ,'LOG-Goldstein and Price'...
    ,'Easom (ES) 2'...
    ,'De Joung (DJ) 3'...
    ,'Rosenbrock 5 (R5) 5'...
    ,'Rosenbrock 10 (R10) 10'...
    ,'Rosenbrock 20 (R20) 20'...
    ,'Rosenbrock 50 (R50) 50'...
    ,'Rosenbrock 100 (R100) 100'...
    ,'Zakharov 2 (Z2) 2'...
    ,'Zakharov 5 (Z5) 5'...
    ,'Zakharov 10 (Z10) 10'...
    ,'Zakharov 20 (Z20) 20'...
    ,'Zakharov 50 (Z50) 50'...
    ,'Zakharov 100 (Z100) 100'...
    ,'Bohachevsky 1 (B1) 2'...
    ,'Bohachevsky 2 (B2) 2'...
    ,'Bohachevsky 3 (B3) 2'...
    ,'Griewank (GR) 6'...
    ,'Aaron 1 2'...
    ,'Aaron 2 2'...
    ,'Booth'...
    ,'Dixon and Price 7'...
    ,'Spiral'...
    ,'EXP-Shekel 5'); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

IntVars = []; uP = []; fGoal = []; nLocal = [];
nGlobal = []; f_Low = []; x_opt = []; f_opt = [];
x_0 = [];
b_L=[]; b_U=[]; A=[];

if P == 1
    Name = 'Shekel 5';
    x_L = [ 0  0  0  0]';
    x_U = [10 10 10 10]';
    x_0 = 0.5*(x_U+x_L);
    x_opt = [4.000036, 4.000132, 4.00037, 4.000132];
    f_opt = -10.1531996784633;
    f_Low=-20;
    x_min = [ 0  0  0  0];
    x_max = [10 10 10 10];
    nGlobal = 1; % Number of global minima
    nLocal  = 5; % Number of local minima
elseif P == 2
    Name = 'Shekel 7';
    x_L = [ 0  0  0  0]';
    x_U = [10 10 10 10]';
    x_0 = []; % Used by EGO for maximizing Likelihood.
    x_opt = [4.000571, 4.000689, 3.999488, 3.999605];
    f_opt = -10.4029405660724;
    f_Low=-20;
    x_min = [ 0  0  0  0];
    x_max = [10 10 10 10];
    nGlobal = 1;
    nLocal  = 7;
elseif P == 3
    Name = 'Shekel 10';
    x_L = [ 0  0  0  0]';
    x_U = [10 10 10 10]';
    x_0 = []; % Used by EGO for maximizing Likelihood.
    x_opt = [4.000747, 4.000593, 3.999663, 3.999510];
    f_opt = -10.5364098166561;
    f_Low=-20;
    x_min = [ 0  0  0  0];
    x_max = [10 10 10 10];
    nGlobal =  1;
    nLocal  = 10;
elseif P == 4
    Name = 'Hartman 3';
    x_L = [0 0 0]';
    x_U = [1 1 1]';
    x_0 = [0.0067 6.012 11.8]'; % Used by EGO for maximizing Likelihood.
    x_opt = [0.114614,0.555648,0.852546];
    f_opt = -3.86278214778846;
    f_Low=-10;
    x_min = [0 0 0];
    x_max = [1 1 1];
    nGlobal = 1;
    nLocal  = 4;
elseif P == 5
    Name = 'Hartman 6';
    x_L = [0 0 0 0 0 0]';
    x_U = [1 1 1 1 1 1]';
    x_0 = [2.6852 0.8333 -0.2778 0.8333 1.0185 -0.0926]'; % Used by EGO for maximizing Likelihood.
    x_opt = [0.201689,0.150010,0.476874,0.275332,0.311651,0.657300];
    f_opt = -3.32236801139709;
    f_Low=-10;
    x_min = [0 0 0 0 0 0];
    x_max = [1 1 1 1 1 1];
    nGlobal = 1;
    nLocal  = 4;
elseif P == 6
    Name = 'Branin RCOS';
    x_L = [-5 0]';
    x_U = [10 15]';
    x_0 = log([0.032  0.0018]); % Used by EGO for maximizing Likelihood.
    x_opt = [ pi 2.275 ; 3*pi 2.475 ; -pi 12.275 ];
    f_opt = 0.397887357729739; % f_opt = 5/(4*pi);
    f_Low=0;
    x_min = [-5 0];
    x_max = [10 15];
    nGlobal = 3;
    nLocal  = 3;
elseif P == 7
    Name = 'Goldstein and Price';
    x_L = [-2 -2]';
    x_U = [ 2  2]';
    x_0 = log([1.9 0.8]'); % Used by EGO for maximizing Likelihood.
    x_opt = [0 -1];
    f_opt = 3;
    f_Low=0;
    x_min = [-2 -2];
    x_max = [ 2  2];
    nGlobal = 1;
    nLocal  = 4;
elseif P == 8
    % Yao Y.
    % "Dynamic Tunneling Algorithm for Global Optimization"
    % IEEE Transactions on systems, Man, and Cybernetics
    % Vol 19, Page 1228, Example 6.
    Name = 'Six-Hump Camel';
    x_L = [-3 -2]';
    x_U = [ 3  2]';
    x_0 = [1.8896 -4.9989 ]'; % Used by EGO for maximizing Likelihood.
    x_opt = [0.089842 -0.712656 ; -0.089842 0.712656];
    f_opt = -1.031628453488552;
    f_Low=-5;
    x_min = [-3 -2];
    x_max = [ 3  2];
    nGlobal = 2;
    nLocal  = 6;
elseif P == 9
    % Yao Y.
    % "Dynamic Tunneling Algorithm for Global Optimization"
    % IEEE Transactions on systems, Man, and Cybernetics
    % Vol 19, Page 1229, Example 7.
    Name = 'Two-Dimensional Shubert';
    x_L = [-10 -10]';
    x_U = [ 10  10]';
    x_0 = []; % Used by EGO for maximizing Likelihood.
    x_opt = [  5.4826   -1.42513
        -0.80032   4.85805
        4.85805  -0.80032
        -7.08350  -7.70831
        -0.80032  -7.70831
        5.48286  -7.70831
        -7.70831  -7.08350
        -1.42513  -7.08350
        4.85805  -7.08350
        -7.08350  -1.42513
        -0.80032  -1.42513
        -7.70831  -0.80032
        -1.42513  -0.80032
        -7.08350   4.85805
        5.48286   4.85805
        -7.70831   5.48286
        -1.42513   5.48286
        4.85805   5.48286  ];
    f_opt = -186.7309088256903;
    f_Low=-200;
    x_min = [-10 -10];
    x_max = [ 10  10];
    nGlobal = 18;
    nLocal  = 760;
elseif P == 10 | P == 11 | P == 12
    if P == 10
        Name = 'Sphere model 2';
        n = 2;
    elseif P == 11
        Name = 'Sphere model 5';
        n = 5;
    elseif P == 12
        Name = 'Sphere model 10';
        n = 10;
    end
    fGoal = 1e-6;
    f_Low = 0;
    f_opt = 0;
    x_opt = ones(1,n);
    x_L = -5*ones(n,1);
    x_U =  5*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 13 | P == 14 | P == 15
    if P == 13
        Name = 'Griewanks function 2';
        n = 2;
    elseif P == 14
        Name = 'Griewanks function 5';
        n = 5;
    elseif P == 15
        Name = 'Griewanks function 10';
        n = 10;
    end
    fGoal = 1e-4;
    f_Low = 0;
    f_opt = 0;
    x_opt = 100*ones(1,n);
    x_L   = -600*ones(n,1);
    x_U   =  600*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 16 | P == 17 | P == 18
    x_opt = [];
    f_opt = [];
    if P == 16
        Name  = 'Shekels foxholes 2';
        n     = 2;
        fGoal = -12.1190083798;
        f_opt = fGoal; % Guess
        x_opt = [8.024065,9.146534]; % from glbSolve / npsol
    elseif P == 17
        Name = 'Shekels foxholes 5';
        n     = 5;
        % from glbSolve / npsol
        fGoal = -10.4039206;
        f_opt = fGoal; % Guess
        x_opt = [8.024917,9.151728,5.113927,7.620861,4.564085];
    elseif P == 18
        %fGoal = -9; 
        Name  = 'Shekels foxholes 10';
        n     = 10;
	% Solution obtained from new version of glcCluster 091205
        %x_opt = [ 8.024725573166975 9.151742042096549 5.113852745456943 7.621090829932100 ...
        %          4.564087033891309 4.710919434871751 2.995996631819100 6.125845669478960 ...
	%	  0.734151702778434 4.982135926165935]';
        %fGoal = -10.207857632486286;
	% SNOPT improved x_opt slightly starting from glcCluster solution
        fGoal = -10.207876836780152;
        x_opt = [ 8.024965290147895 9.151926376025360 5.113988849847106 7.620956984993030 4.564018471706530 ...
                  4.711003025383680 2.996029082447810 6.125990970569745 0.734056774552124 4.981997675469230];
	
        f_opt = fGoal; % Guess
	
	    
    end
    f_Low = -20;
    x_L   = zeros(n,1);
    x_U   = 10*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 19 | P == 20 | P == 21
    if P == 19
        Name  = 'Michalewiczs function 2';
        fGoal = -1.80130341008983;
        f_opt = fGoal; % Guess
        x_opt = [2.202905,1.570796]; % from glbSolve / npsol
        n     = 2;
    elseif P == 20
        Name  = 'Michalewiczs function 5';
        fGoal = -4.687658179;  % from glbSolve + npsol
        f_opt = fGoal; % Guess
        x_opt = [2.202906,1.570796,1.284992,1.923058,1.720470];
        n     = 5;
    elseif P == 21
        Name  = 'Michalewiczs function 10';
        fGoal = -7.872489176499680; % from glcCluster with NPSOL
        fGoal = -8.111637967813923; % from multiMin with SNOPT
        x_opt = [ 0.073258687813626 1.570796281349359 2.219333353327856 1.113781415521201 1.720469814075476 ...
                  1.570796403991472 2.221053309845177 1.756086528146987 1.958964957372832 2.329651907056382 ];
        n     = 10;
        f_opt = fGoal; % Guess
    end
    x_L   = zeros(n,1);
    x_U   = pi*ones(n,1);
    f_Low = -20;
    x_min = x_L;
    x_max = x_U;
elseif P == 22 | P == 23 | P == 24
    f_opt = [];
    x_opt = [];
    if P == 22
        Name  = 'Langermans function 2';
        n     = 2;
        fGoal = -3.0677475617235372; % Accuracy from snopt
        f_opt = fGoal; % Guess
        x_opt = [7.735257,8.860257];
    elseif P == 23
        Name  = 'Langermans function 5';
        n     = 5;
        fGoal = -1.499943823715334600; % From snopt
        f_opt = fGoal; % Guess
        x_opt = [ 8.02329762255323 9.15653111755802 5.11632952767342 ...
            7.61983157672815 4.56378371087901];
    elseif P == 24
        Name  = 'Langermans function 10';
        n     = 10;
        fGoal = -0.96399998152060;
        f_opt = fGoal; % Guess
        x_opt = [ 6.30598735426646 8.58299141102869 6.08399254062668 ...
            1.13800014809038 4.35000836618424 3.13399328043645 ...
            7.85300210477163 6.06100677670264 7.45700606167984 ...
            2.25799751564348 ];
    end
    f_Low = -5;
    x_L   = zeros(n,1);
    x_U   = 10*ones(n,1);
    x_min = x_L;
    x_max = x_U;
elseif P == 25
    Name  = 'HGO 468:1';
    x_L   = [0 0]';
    x_U   = [1 1]';
    x_0   = [];
    x_opt = [1 0.63492204];
    f_opt = -2.51997258;
    f_Low = -5;
    x_min = [0 0];
    x_max = [1 1];
elseif P == 26
    Name  = 'HGO 468:2';
    x_L   = [0 0]';
    x_U   = [1 1]';
    x_0   = [];
    x_opt = [0.28539815 0];
    f_opt = -2.81859485;
    f_Low = -5;
    x_min = [0 0];
    x_max = [1 1];
elseif P == 27
    Name  = 'Hock-Schittkowski 5';
    x_L   = [-1.5 -3]';
    x_U   = [4 3]';
    x_opt = [-pi/3+0.5,-pi/3-0.5];
    f_opt = -0.5*sqrt(3)-pi/3;
    f_Low = -5;
    x_min = x_L;
    x_max = x_U;
    x_0 = [0 0]'; % If running local solver
    nGlobal = [];
    nLocal  = [];
elseif P == 28
    Name  = 'RB BANANA';
    x_L   = [-2 -2]';
    x_U   = [ 2  2]';
    x_0   = [-1.9604 -4.9996]'; % Used by EGO for maximizing Likelihood.
    x_opt = [1 1];
    f_opt = 1;
    f_Low = 0 ;
    x_min = [-2 -2];
    x_max = [ 2  2];
    nGlobal = 1;
    nLocal  = 1;
elseif P == 29	% Myers smooth
    Name = 'Myers smooth';
    x_L = [0.5 0.5]';
    x_U = [3.5 3.5]';
    x_opt = [3.5 3.5];
    f_opt = 5.715815695448042;
    f_Low = 0;
    x_min = x_L;
    x_max = x_U;
    nGlobal = [];
    nLocal  = [];
    x_0 = [2 2]';
elseif P == 30	% Myers smoothly fluctuating
    Name = 'Myers smoothly fluctuating';
    x_L = [0.5 0.5]';
    x_U = [3.5 3.5]';
    x_opt = [0.500229 2.510288 ;
        2.510288 0.500229];
    f_opt = 5.8861710434083498;
    f_Low = 0;
    x_min = x_L;
    x_max = x_U;
    nGlobal = [];
    nLocal  = [];
    x_0 = [2 2]';
elseif P == 31	% Myers strongly fluctuating
    Name = 'Myers strongly fluctuating';
    x_L = [0.5 0.5]';
    x_U = [3.5 3.5]';
    x_opt = [3.358939 0.500076 ;
        0.500076 3.358939];
    f_opt = 5.508607992798803;
    f_Low = 0;
    x_min = x_L;
    x_max = x_U;
    nGlobal = [];
    nLocal  = [];
    x_0 = [2 2]';
elseif P == 32
    Name  = 'LOG-Goldstein and Price';
    x_L   = [-2 -2]';
    x_U   = [ 2  2]';
    x_0   = log([1.9 0.8]'); % Used by EGO for maximizing Likelihood.
    x_opt = [0 -1];
    f_opt = log(3);
    f_Low = 0;
    x_min = [-2 -2];
    x_max = [ 2  2];
    nGlobal = 1;
    nLocal  = 4;
elseif P == 33	% Easom (ES) 2
    Name  = 'Easom (ES) 2';
    x_L   = [-10 -10]';
    x_U   = [10 10]';
    % Bounds 10 according to M.M. Ali et al. JOGO (2005) 31:635-672
    %x_L   = [-100 -100]';
    %x_U   = [100 100]';
    x_opt = [pi pi];
    f_opt = -1;
    x_min = x_L;
    x_max = x_U;
elseif P == 34	% De Joung (DJ) 3
    Name = 'De Joung (DJ) 3';
    x_L = [-5.12 -5.12 -5.12]';
    x_U = [5.12 5.12 5.12]';
    x_opt = [0 0 0];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 35 | P == 36 | P == 37 | P == 38 | P == 39
    if P == 35
        Name = 'Rosenbrock 5 (R5) 5';
        n = 5;
    elseif P == 36
        Name = 'Rosenbrock 10 (R10) 10';
        n = 10;
    elseif P == 37
        Name = 'Rosenbrock 20 (R20) 20';
        n = 20;
    elseif P == 38
        Name = 'Rosenbrock 50 (R50) 50';
        n = 50;
    elseif P == 39
        Name = 'Rosenbrock 100 (R100) 100';
        n = 100;
    end
    x_L = -5*ones(n,1);
    x_U = 10*ones(n,1);
    x_opt = ones(1,n);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 40 | P == 41 | P == 42 | P == 43 | P == 44 | P == 45
    if P == 40
        Name = 'Zakharov 2 (Z2) 2';
        n = 2;
    elseif P == 41
        Name = 'Zakharov 5 (Z5) 5';
        n = 5;
    elseif P == 42
        Name = 'Zakharov 10 (Z10) 10';
        n = 10;
    elseif P == 43
        Name = 'Zakharov 20 (Z20) 20';
        n = 20;
    elseif P == 44
        Name = 'Zakharov 50 (Z50) 50';
        n = 50;
    elseif P == 45
        Name = 'Zakharov 100 (Z100) 100';
        n = 100;
    end
    x_L = -5*ones(n,1);
    x_U = 10*ones(n,1);
    x_opt = zeros(1,n);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 46	% Bohachevsky 1 (B1) 2
    % M.M. Ali, C. Khompatraporn and Z.B. Zabinsky, "A numerical
    % evaluation of several stochastic algorithms on selected continuous
    % global optimization test problems", Journal of Global Optimization
    % (2005) 31:635-672
    Name = 'Bohachevsky 1 (B1) 2';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_opt = [0 0];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 47	% Bohachevsky 2 (B2) 2
    Name = 'Bohachevsky 2 (B2) 2';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_opt = [0 0];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 48	% Bohachevsky 3 (B3) 2
    Name = 'Bohachevsky 3 (B3) 2';
    x_L = [-10 -10]';
    x_U = [10 10]';
    f_opt = 0;
    x_opt = [0 0];
    x_min = x_L;
    x_max = x_U;
elseif P == 49	% Griewank (GR) 6
    Name = 'Griewank (GR) 6';
    x_L = -1*ones(6,1);
    x_U = -x_L;
    x_opt = zeros(1,6);
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 50	% Aaron 1 2
    Name = 'Aaron 1 2';
    x_L = [-1 -1]';
    x_U = -x_L;
    x_opt = [0 0];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 51	% Aaron 2 2
    Name = 'Aaron 2 2';
    x_L = [-1 -1]';
    x_U = -x_L;
    x_opt = [-1 -1 ;
        1 -1];
    f_opt = -5.913959886547349;
    x_min = x_L;
    x_max = x_U;
elseif P == 52	% Booth
    Name = 'Booth';
    x_L   = [-10;-10];
    x_U   = [ 10; 10];
%     x_opt = [131,32]/36;
%     f_opt = 1805/144;
    x_opt = [1,3];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 53	% Dixon & Price
    Name = 'Dixon and Price';
    x_L   = [-10;-10];
    x_U   = [ 10; 10];
    x_opt = [1 1/sqrt(2) ; 1 -1/sqrt(2)];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 54
    Name='Spiral function';
    uP(2,1) = 0.4;
    uP(1,1) = 1;
    x_L   = [-10;-10];
    x_U   = [ 10; 10];
    x_0   = []; % Used by EGO for maximizing Likelihood.
    x_opt = [0,-1E-10];     % x_opt = [0,0-], ie. lim x(2) approaches 0 from the left.
    f_opt = -1;
    f_Low = - 100;
    x_max = [ 10  10];
    x_min = [-10 -10];
elseif P == 55
    Name = 'EXP-Shekel 5';
    x_L = [ 0  0  0  0]';
    x_U = [10 10 10 10]';
    x_0 = 0.5*(x_U+x_L);
    x_opt = [4.000036, 4.000132, 4.00037, 4.000132];
    f_opt = 1E5*exp(-10.1531996784633);
    f_Low = 0;
    x_min = [ 0  0  0  0];
    x_max = [10 10 10 10];
    nGlobal = 1;
    nLocal  = 5;
else
    error('glb_prob: Illegal problem number');
end
% Set x_0 to zeros (dummy for GUI)
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the Prob
Prob = glcAssign('glb_f', x_L, x_U, Name, A, b_L, b_U, ...
    [], [], [], x_0, IntVars, [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.uP = uP;
Prob.MIP.nLocal  = nLocal;
Prob.MIP.nGlobal = nGlobal;
Prob.MIP.fGoal=fGoal(:);

% MODIFICATION LOG:
%
% 980915  mbk  Known optimal values added.
% 980922  hkh  Change name f_min to f_Low
% 980923  mbk  Check if x_0 is empty, if so then set x_0 to origo.
% 981005  mbk  b_L=[]; b_U=[]; c_L=[]; c_U=[]; for all problems.
% 981006  hkh  Added call to checkuP
% 981011  hkh  Changed to use tomFiles for name definitions
% 981022  hkh  Set Name=[] if P is illegal
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 981127  mbk  x_0 defined for some problems, used by EGO for the
%              maximization of the Likelihood.
% 990413  mbk  Added problem Hock-Schittkowski 5
% 990623  hkh  Use glcVarDef and glcProbSet instead of ProbVarDef and ProbSet
% 041117  med  xxx_prob removed and code added
% 050405  med  Added 19 problems
% 050423  hkh  Added x_opt to Shekel problems
% 080420  hkh  Add suggested f_opt=fGoal, switch LOG-GP--spiral, 1-51 runnable
% 080423  nhq  Added new x_opt values and some alternative x_opt points
% 080423  nhq  More decimals in some optimal f(x)
% 080603  med  Switched to glcAssign, cleaned
% 080915  nhq  Added x_opt and f_opt to problems 33,52,53,54.
% 080918  nhq  Improved x_opt and f_opt for problem 8.
% 080925  nhq  Improved x_opt and f_opt for problem 52.
% 091205  hkh  Improved x_opt and f_opt for problem 18 and 21

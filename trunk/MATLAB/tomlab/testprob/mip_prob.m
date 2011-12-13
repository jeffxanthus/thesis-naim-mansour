% mip_prob:
%
% Defines Mixed-Integer Programming problems
%
% function [probList, Prob] = mip_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc. $Release: 7.3.0$
% Written June 1, 1999.   Last modified Oct 15, 2008.

function [probList, Prob] = mip_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    '1st simple LP. OR course'...
    ,'Luenberger 2.6-2'...
    ,'Example LP-3 OA 5p'...
    ,'Luenberger 2.6-3'...
    ,'Luenberger 2.6-5'...
    ,'Marshall-Suurballe cycling example'...
    ,'Weingartner 1 - 2/28 0-1knapsack'...
    ,'Kuhn example.'...
    ,'Beale cycling example.'...
    ,'Winston Ch.6 Review 17. Dual feasible'...
    ,'Fletcher 8.21.3.9. Redundancy in constr.'...
    ,'Winston Ex. 4.12 B4. Max | |. Rewritten'...
    ,'Hansen, Plateau 1 - 4/28 0-1 knapsack'...
    ,'PB 4 - 2/29 0-1 knapsack'...
    ,'Generalized Assignment, Wolsey-98 pp159'...
    ,'Generalized Assignment, Wolsey-98 9.8.16'...
    ,'University of Washington thrust problem'...
    ,'aflow30a'...
    ,'bienst1'...
    ,'bienst2'...
    ,'danoint'...
    ,'fixnet6'...
    ,'glass4'...
    ,'markshare1'...
    ,'markshare1_1'...
    ,'markshare2'...
    ,'markshare2_1'...
    ,'mas74'...
    ,'mas76'...
    ,'modglob'...
    ,'neos14'...
    ,'neos15'...
    ,'neos16'...
    ,'neos5'...
    ,'noswot'...
    ,'opt1217'...
    ,'pk1'...
    ,'pp08a'...
    ,'pp08aCUTS'...
    ,'qiu'...
    ,'ran14x18_1'...
    ,'rout'...
    ,'set1ch'...
    ,'timtab1'...
    ,'timtab2'...
    ,'tr12-30'...
    ,'vpm2'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_0 = []; VarWeight = []; KNAPSACK = []; x_opt = []; f_opt = [];
optParam = optParamDef; f_Low = []; x_min = []; x_max = [];
Alg = 0; c = []; A = []; b_U = [];

if P == 1
    Name='1st simple LP. OR course';
    % Problem formulated on standard form
    A   = [1 2 1 0
        4 1 0 1 ];
    b_U = [6;12];
    c   = [-7  -5  0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = [2 1 6 12]';
    x_opt=[2 1 2 3];
    f_opt=-19;
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=3*ones(n,1);
    IntVars = 1:4;
elseif P == 2
    Name='Luenberger 2.6-2';
    % Problem formulated on standard form
    % Shows that it is unnessary to add sum(x(i)=1,
    % because summing the two first constraints gives sum(x(i))=1
    A=[0.1 0.25 0.5 0.75 0.95; 0.9 0.75 0.5 0.25 0.05;ones(1,5) ];
    b_U=[0.3 0.7 1]';
    c   = [5  4 3 2 1.5]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 3
    Name='Example LP-3 from OA 5p';
    % Problem formulated on standard form
    A1=[ones(1,3),zeros(1,3);zeros(1,3),ones(1,3)];
    A2=[eye(3),eye(3);[6 0.4 0.26],zeros(1,3);zeros(1,3),[12/5 4/25 13/125]];
    A3=[12 7.5 10 0 0 0;0 0 0 24/5 3 4];
    A  = [A1,zeros(2,7);A2,eye(5),zeros(5,2);A3,zeros(2,5),-eye(2)];
    b_U= [10 25 12 15 25 10 8 90 95]';
    c  = [1500 2400 3000 1500 2400 3000,zeros(1,7)]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    % If running DIRECT
    x_U(1:4) = 5;
    x_U(5:n) = 20;
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = (1:n-size(A,1));
elseif P == 4
    Name='Luenberger 2.6-3';
    % Problem formulated on standard form
    A  = [0.3 0 0 0.3 0 0;0 0.2 0 0 0.4 0; 0 0 0.3 0 0 0.2];
    b_U= [900000 800000 500000 ]';
    A=A/1E0;
    b_U=b_U/1E0;
    c  = [35 35 35 30 30 30]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 5
    Name='Luenberger 2.6-5';
    % Problem formulated on standard form, by adding two vars
    A   = [2 -2 -2 1; 1 -1 0 -1];
    b_U = [4 1]';
    c   = [1 -1 4 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
    optParam.MaxIter=5000;
elseif P == 6
    Name='Marshall-Suurballe cycling example.';
    A   = -[-.5 5.5 2.5 -9  -1 0 0
        -.5 1.5 0.5 -1  0 -1 0
        -1 0   0    0  0 0 -1];
    b_U = [0 0 1]';
    c   = -[10 -57 -9 -24 0 0 0]';
    %B=~[0 0 0 0 1 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 7
    Name='Weingartner 1 - 2/28 0-1 knapsack';
    % Problem formulated as a minimum problem
    A = [ 45      0     85     150     65     95     30      0    170  0 ...
        40     25     20       0      0     25      0      0     25  0 ...
        165      0     85       0      0      0      0    100  ; ...
        30     20    125       5     80     25     35     73     12  15 ...
        15     40      5      10     10     12     10      9      0  20 ...
        60     40     50      36     49     40     19    150];
    b_U = [600;600];  % // 2 knapsack capacities
    c   = [1898  440  22507  270  14148   3100   4650  30800   615   4975 ...
        1160   4225    510   11880    479    440    490    330    110    560 ...
        24355   2885  11748    4550    750   3720   1950  10500]';  % 28 weights
    [Prob,A,b_U,c,x_0,x_L,x_U,n,IntVars,...
        KNAPSACK, VarWeight, b_L,x_min,x_max]=knapmake(ProbDef,A,b_U,c);
    Alg = 2;
    optParam.MaxIter = 5000;
    f_opt=-141278;
    f_Low=-1E7;
    % References for 0/1 Multiple Knapsack Problems, solved in:
    % Tabu Search:
    % F. Dammeyer and S. Voss (1991) "Dynamic tabu list management using
    % the reverse elimination method."
    % Annals of Operations Research.

    % Simulated Annealing:
    % A. Drexel (1988) "A Simulated Annealing Approach to the Multiconstraint
    % Zero-One Knapsack Problem."
    % Computing, 40:1-8.

    % Genetic Algorithms:
    % S. Khuri, T. Baeck, J. Heitkoetter, (1994) "The Zero/One Multiple
    % Knapsack Problem and Genetic Algorithms", to appear in the 1994
    % ACM Symposium on Applied Computing, SAC'94, Phoenix, Arizona.
elseif P == 8
    Name='Kuhn example.';
    A   = [-2 -9 1 9 1 0
        1/3 1 -1/3 -2 0 1];
    b_U = [0 0 ]';
    c   = [-2 -3 1 12  0 0]';
    %B   = ~[0 0 0 0 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 9
    % Bazaraa page 165
    Name='Beales cycling example.';
    A = [eye(3), [0.25 -8 -1 9; 0.5 -12 -0.5 3; 0 0 1 0]];
    b_U = [0 0 1]';
    c   = [0 0 0 -0.75 20 -0.5 6]';
    %B   = ~[1 1 1 0 0 0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 10
    Name='Winston Ch.6 Review 17. Dual feasible';
    % Could be solved by dual LP
    A   = [1  1 -1  0
        1 -2  0 -1 ];
    b_U = [5 8]';
    c   = [2 1 0 0]';
    %B= ~[1 1 0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 11
    Name='Fletcher 8.21. Redundancy in constr.';
    A   = [1  0  2 1 0
        0  1 -1 0 1
        1  1  1 0 0 ];
    b_U = [1  1 2 ]';
    c   = [-1 -1 0 0 0]';
    if 0
        A=A(1:2,:);
        b_U=b_U(1:2);
    end
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:n;
elseif P == 12
    Name='Winston Ex. 4.12 B4. Max | |. Rewritten';
    A   = [4    1    0 0 1 0
        2   -1    0 0 0 1
        2   -3   -1 1 0 0
        ];
    b_U = [4  0.5 0]';
    c   = [0 0  -1 -1 0 0]';
    %B   =~ [0 0 0 1 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
    IntVars = 1:3;
elseif P == 13
    Name='Hansen, Plateau 1 - 4/28 0-1 knapsack';
    % Problem formulated as a minimum problem
    A = [ 40	 91	3    12     3	 18    25     1     1	  8 ...
        1	  1    49     8    21	  6	1     5     8	  1 ...
        0	 42	6     4     8	  0    10     1 ;     ...
        16	 92	4    18     6	  0	8     2     1	  6 ...
        2	  1    70     9    22	  4	1     5     6	  0 ...
        4	  8	4     3     0	 10	0     6 ;     ...
        38	 39	5    40     8	 12    15     0     1	 20 ...
        3	  0    40     6     8	  0	6     4     4	  1 ...
        5	  8	2     8     0	 20	0     0 ;     ...
        38	 52	7    20     0	  3	4     1     2	  4 ...
        6	  1    18    15    38	 10	4     8     0	  3 ...
        0	  6	1     3     0	  3	5     4 ];
    b_U = [219	203   208   180]';
    c   = [ 560  1125    68   328    47	122   196    41    25	115 ...
        82	 22   631   132   420	 86    42   103    81	 26 ...
        49	316    72    71    49	108   116    90]';
    [Prob,A,b_U,c,x_0,x_L,x_U,n,IntVars,...
        KNAPSACK, VarWeight, b_L,x_min,x_max]=knapmake(ProbDef,A,b_U,c);
    Alg = 2;
    optParam.MaxIter = 5000;
    f_opt=-3418;
    f_Low=-1E7;
elseif P == 14
    Name='PB 4 - 2/29 0-1 knapsack';
    % Problem formulated as a minimum problem
    A = [ 25  17  20  22 15 10  50 10 150  0  0 0 40 60 , zeros(1,15); ...
        zeros(1,5), 2 5 6 40 2	6 10 13 30 15 5 5 10 15	91 24 15 ...
        15 5 10 15 10 10 10];
    b_U = [ 153	 154]';
    c   = [ 7074	5587  5500   3450   367  4235	9262  6155  32240   1600 ...
        4500	6570  7010  16020  8000  2568	2365  4350   4975  29400 ...
        7471	3840  3575   1125  1790  2450	 540   445    220]';
    [Prob,A,b_U,c,x_0,x_L,x_U,n,IntVars, ...
        KNAPSACK, VarWeight, b_L,x_min,x_max]=knapmake(ProbDef,A,b_U,c);
    Alg = 2;
    optParam.MaxIter = 5000;
    f_opt=-95168;
    f_Low=-1E7;
elseif P == 15
    % A generalized assignment problem from Wolsey 1998,
    % Section 9.6, pp159.
    %
    % Define the linear sos1 constraints explicitly
    Name='Generalized Assignment, Wolsey-98 pp159';
    A=[95  1 21 66 59; ...
        54 53 44 26 60; ...
        3 91 43 42  5; ...
        72 30 56 72  9; ...
        44  1 71 13 27; ...
        20 99 87 52 85; ...
        72 96 97 73 49; ...
        75 82 83 44 59; ...
        68  8 87 74  4; ...
        69 83 98 88 45];
    C=-[110  16 25 78 59; ...
        65  69 54 28 71; ...
        19  93 45 45  9; ...
        89  31 72 83 20; ...
        62  17 77 18 39; ...
        37 115 87 59 97; ...
        89 102 98 74 61; ...
        78  96 87 55 77; ...
        74  27 99 91  5; ...
        88  97 99 99 51];
    b = [91 87 109 88 64]';
    m   = length(b);
    [c,x_L,x_U,b_L,b_U,A]=abc2gap(A,b,C,0);
%if 0                                      % Not necessary to add slacks
%    % Add slacks. First m are inequalities
%    mA       = size(A,1);
%    A        = [A,[eye(m);zeros(mA-m,m)]];
%    x_L      = [x_L;zeros(m,1)];
%    x_U      = [x_U;b];
%    c        = [c; zeros(m,1)];
%    % b_L(1:5) = b_U(1:5);
%end
    IntVars  = (1:length(c));
    f_Low    = -1000;                      % f_Low <= f_optimal must hold
    f_opt    = -545;
    % Set priority weights on the variables
    VarWeight = c;
    optParam.MaxIter = 25000;
    % Tell mipSolve to use knapsack heuristic
    KNAPSACK = 0;

elseif P == 16
    % A generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
    %
    % Define the linear sos1 constraints explicitly
    Name='Generalized Assignment, Wolsey-98 9.8.16';
    A=[ 3 24 53 27 17; ...
        15 23 43 74 23; ...
        54 43 27 21 36; ...
        92 83 45 35 23; ...
        19 10 33 43 12; ...
        91 55 32 26 23; ...
        15 25 35 37 28; ...
        47 43 35 32 37; ...
        34 23 52 46 43; ...
        35 23 34 25 40];
    C=-[15 44 76 43 34; ...
        19 23 45 46 34; ...
        10  6  3 23 15; ...
        60 45 34 36 23; ...
        21 12 34 44 10; ...
        67 35 34 20 37; ...
        23 34 44 47 32; ...
        23 25 32 15 27; ...
        15 13 23 24 34; ...
        10 15 23 12 13];
    b = [80 63 75 98 59]';
    m   = length(b);
    [c,x_L,x_U,b_L,b_U,A]=abc2gap(A,b,C,0);
%if 0                                      % Not necessary to add slacks
%    % Add slacks. First m are inequalities
%    mA  = size(A,1);
%    A   = [A,[eye(m);zeros(mA-m,m)]];
%    x_L = [x_L;zeros(m,1)];
%    x_U = [x_U;b];
%    c   = [c; zeros(m,1)];
%end
    IntVars=(1:length(c));
    f_Low=-1000;              % f_Low <= f_optimal must hold
    optParam.MaxIter=25000;
    % Set priority weights on the variables
    VarWeight=c;
    % Tell mipSolve to use knapsack heuristic
    KNAPSACK=0;
elseif P == 17
    Name='Thrust-problem, University of Washington';
    load('tp');
    x_0=[];
    [m,n] = size(A);
    x_L   = zeros(n,1);
    x_U   = ones(n,1);
    ISlack = eye(m);
    % Change sign on equations with rhs < 0
    ix      = find(b_U < 0);
    A(ix,:) = -A(ix,:);
    b_U(ix) = -b_U(ix);
    ISlack(ix,:) = -ISlack(ix,:);
    % Add slacks
    A   = [A,ISlack];
    x_L = [x_L;zeros(m,1)];
    x_U = [x_U;b_U];
    c   = [c; zeros(m,1)];
    x_min = x_L;
    x_max = x_U;
    b_L = b_U;
    % All original variables should be integer
    IntVars=1:n;
    optParam.MaxIter=25000;
    % Depth First, then Breadth
    Alg=2;
    f_Low=0;        % f_Low <= f_optimal must hold
elseif P == 18
    load mip_probmat P1;
    Name = P1.Name;
    c    = P1.c;
    A    = P1.A;
    b_L  = P1.b_L;
    b_U  = P1.b_U;
    x_L  = P1.x_L;
    x_U  = P1.x_U;
    IntVars = P1.IntVars;
elseif P == 19
    load mip_probmat P2;
    Name = P2.Name;
    c    = P2.c;
    A    = P2.A;
    b_L  = P2.b_L;
    b_U  = P2.b_U;
    x_L  = P2.x_L;
    x_U  = P2.x_U;
    IntVars = P2.IntVars;
elseif P == 20
    load mip_probmat P3;
    Name = P3.Name;
    c    = P3.c;
    A    = P3.A;
    b_L  = P3.b_L;
    b_U  = P3.b_U;
    x_L  = P3.x_L;
    x_U  = P3.x_U;
    IntVars = P3.IntVars;
elseif P == 21
    load mip_probmat P4;
    Name = P4.Name;
    c    = P4.c;
    A    = P4.A;
    b_L  = P4.b_L;
    b_U  = P4.b_U;
    x_L  = P4.x_L;
    x_U  = P4.x_U;
    IntVars = P4.IntVars;
elseif P == 22
    load mip_probmat P5;
    Name = P5.Name;
    c    = P5.c;
    A    = P5.A;
    b_L  = P5.b_L;
    b_U  = P5.b_U;
    x_L  = P5.x_L;
    x_U  = P5.x_U;
    IntVars = P5.IntVars;
elseif P == 23
    load mip_probmat P6;
    Name = P6.Name;
    c    = P6.c;
    A    = P6.A;
    b_L  = P6.b_L;
    b_U  = P6.b_U;
    x_L  = P6.x_L;
    x_U  = P6.x_U;
    IntVars = P6.IntVars;
elseif P == 24
    load mip_probmat P7;
    Name = P7.Name;
    c    = P7.c;
    A    = P7.A;
    b_L  = P7.b_L;
    b_U  = P7.b_U;
    x_L  = P7.x_L;
    x_U  = P7.x_U;
    IntVars = P7.IntVars;
elseif P == 25
    load mip_probmat P8;
    Name = P8.Name;
    c    = P8.c;
    A    = P8.A;
    b_L  = P8.b_L;
    b_U  = P8.b_U;
    x_L  = P8.x_L;
    x_U  = P8.x_U;
    IntVars = P8.IntVars;
elseif P == 26
    load mip_probmat P9;
    Name = P9.Name;
    c    = P9.c;
    A    = P9.A;
    b_L  = P9.b_L;
    b_U  = P9.b_U;
    x_L  = P9.x_L;
    x_U  = P9.x_U;
    IntVars = P9.IntVars;
elseif P == 27
    load mip_probmat P10;
    Name = P10.Name;
    c    = P10.c;
    A    = P10.A;
    b_L  = P10.b_L;
    b_U  = P10.b_U;
    x_L  = P10.x_L;
    x_U  = P10.x_U;
    IntVars = P10.IntVars;
elseif P == 28
    load mip_probmat P11;
    Name = P11.Name;
    c    = P11.c;
    A    = P11.A;
    b_L  = P11.b_L;
    b_U  = P11.b_U;
    x_L  = P11.x_L;
    x_U  = P11.x_U;
    IntVars = P11.IntVars;
elseif P == 29
    load mip_probmat P12;
    Name = P12.Name;
    c    = P12.c;
    A    = P12.A;
    b_L  = P12.b_L;
    b_U  = P12.b_U;
    x_L  = P12.x_L;
    x_U  = P12.x_U;
    IntVars = P12.IntVars;
elseif P == 30
    load mip_probmat P13;
    Name = P13.Name;
    c    = P13.c;
    A    = P13.A;
    b_L  = P13.b_L;
    b_U  = P13.b_U;
    x_L  = P13.x_L;
    x_U  = P13.x_U;
    IntVars = P13.IntVars;
elseif P == 31
    load mip_probmat P14;
    Name = P14.Name;
    c    = P14.c;
    A    = P14.A;
    b_L  = P14.b_L;
    b_U  = P14.b_U;
    x_L  = P14.x_L;
    x_U  = P14.x_U;
    IntVars = P14.IntVars;
elseif P == 32
    load mip_probmat P15;
    Name = P15.Name;
    c    = P15.c;
    A    = P15.A;
    b_L  = P15.b_L;
    b_U  = P15.b_U;
    x_L  = P15.x_L;
    x_U  = P15.x_U;
    IntVars = P15.IntVars;
elseif P == 33
    load mip_probmat P16;
    Name = P16.Name;
    c    = P16.c;
    A    = P16.A;
    b_L  = P16.b_L;
    b_U  = P16.b_U;
    x_L  = P16.x_L;
    x_U  = P16.x_U;
    IntVars = P16.IntVars;
elseif P == 34
    load mip_probmat P17;
    Name = P17.Name;
    c    = P17.c;
    A    = P17.A;
    b_L  = P17.b_L;
    b_U  = P17.b_U;
    x_L  = P17.x_L;
    x_U  = P17.x_U;
    IntVars = P17.IntVars;
elseif P == 35
    load mip_probmat P18;
    Name = P18.Name;
    c    = P18.c;
    A    = P18.A;
    b_L  = P18.b_L;
    b_U  = P18.b_U;
    x_L  = P18.x_L;
    x_U  = P18.x_U;
    IntVars = P18.IntVars;
elseif P == 36
    load mip_probmat P19;
    Name = P19.Name;
    c    = P19.c;
    A    = P19.A;
    b_L  = P19.b_L;
    b_U  = P19.b_U;
    x_L  = P19.x_L;
    x_U  = P19.x_U;
    IntVars = P19.IntVars;
elseif P == 37
    load mip_probmat P20;
    Name = P20.Name;
    c    = P20.c;
    A    = P20.A;
    b_L  = P20.b_L;
    b_U  = P20.b_U;
    x_L  = P20.x_L;
    x_U  = P20.x_U;
    IntVars = P20.IntVars;
elseif P == 38
    load mip_probmat P21;
    Name = P21.Name;
    c    = P21.c;
    A    = P21.A;
    b_L  = P21.b_L;
    b_U  = P21.b_U;
    x_L  = P21.x_L;
    x_U  = P21.x_U;
    IntVars = P21.IntVars;
elseif P == 39
    load mip_probmat P22;
    Name = P22.Name;
    c    = P22.c;
    A    = P22.A;
    b_L  = P22.b_L;
    b_U  = P22.b_U;
    x_L  = P22.x_L;
    x_U  = P22.x_U;
    IntVars = P22.IntVars;
elseif P == 40
    load mip_probmat P23;
    Name = P23.Name;
    c    = P23.c;
    A    = P23.A;
    b_L  = P23.b_L;
    b_U  = P23.b_U;
    x_L  = P23.x_L;
    x_U  = P23.x_U;
    IntVars = P23.IntVars;
elseif P == 41
    load mip_probmat P24;
    Name = P24.Name;
    c    = P24.c;
    A    = P24.A;
    b_L  = P24.b_L;
    b_U  = P24.b_U;
    x_L  = P24.x_L;
    x_U  = P24.x_U;
    IntVars = P24.IntVars;
elseif P == 42
    load mip_probmat P25;
    Name = P25.Name;
    c    = P25.c;
    A    = P25.A;
    b_L  = P25.b_L;
    b_U  = P25.b_U;
    x_L  = P25.x_L;
    x_U  = P25.x_U;
    IntVars = P25.IntVars;
elseif P == 43
    load mip_probmat P26;
    Name = P26.Name;
    c    = P26.c;
    A    = P26.A;
    b_L  = P26.b_L;
    b_U  = P26.b_U;
    x_L  = P26.x_L;
    x_U  = P26.x_U;
    IntVars = P26.IntVars;
elseif P == 44
    load mip_probmat P27;
    Name = P27.Name;
    c    = P27.c;
    A    = P27.A;
    b_L  = P27.b_L;
    b_U  = P27.b_U;
    x_L  = P27.x_L;
    x_U  = P27.x_U;
    IntVars = P27.IntVars;
elseif P == 45
    load mip_probmat P28;
    Name = P28.Name;
    c    = P28.c;
    A    = P28.A;
    b_L  = P28.b_L;
    b_U  = P28.b_U;
    x_L  = P28.x_L;
    x_U  = P28.x_U;
    IntVars = P28.IntVars;
elseif P == 46
    load mip_probmat P29;
    Name = P29.Name;
    c    = P29.c;
    A    = P29.A;
    b_L  = P29.b_L;
    b_U  = P29.b_U;
    x_L  = P29.x_L;
    x_U  = P29.x_U;
    IntVars = P29.IntVars;
elseif P == 47
    load mip_probmat P30;
    Name = P30.Name;
    c    = P30.c;
    A    = P30.A;
    b_L  = P30.b_L;
    b_U  = P30.b_U;
    x_L  = P30.x_L;
    x_U  = P30.x_U;
    IntVars = P30.IntVars;
else
    error('mip_prob: Illegal problem number');
end

Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
    [], [], IntVars, VarWeight, KNAPSACK, [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);

Prob.P          = P;
Prob.Solver.Alg = Alg;
Prob.optParam   = optParam;

function [Prob, AA, bb, cc, x_0, x_L, x_U, n, IntVars, ...
    KNAPSACK, VarWeight, b_L, x_min, x_max] ...
    = knapmake(Prob, A, b, c)

% Knapsack demonstration examples for Integer Programming
%
% knapmake defines the problem in the TOMLAB format
% Make problem on standard form for mipSolve, if changing if 0 to if 1
% (Not needed for mipSolve)

[m,n]        = size(A);
if 0
   AA        = [A eye(m,m)];
   % Change sign to make minimum problem
   cc        = [-c;zeros(m,1)];
   bb        = b;
   [mm nn]   = size(AA);
   x_L       = zeros(nn,1);
   x_U       = [ones(n,1);b];
   x_0       = [zeros(n,1);b];
   % All (original) variables should be integer
   IntVars   = 1:n;
   n         = length(cc);
   b_L       = b;
else
   AA        = A;
   % Change sign to make minimum problem
   cc        = -c;
   x_L       = zeros(n,1);
   x_U       = ones(n,1);
   x_0       = zeros(n,1);
   bb        = b;
   b_L       = -inf*ones(m,1);
   IntVars   = 1:n;
end


% Set priority weights on the variables
VarWeight = cc;

% Tell mipSolve to use knapsack heuristic
KNAPSACK  = 1;
x_min     = x_L;
x_max     = x_U;

% MODIFICATION LOG:
%
% 990929 hkh  Added three knapsack problems
% 011030 hkh  Changed problem names, avoiding too long names
% 031125 ango Fixed x_opt, f_opt for Problem 1
% 041110 med  Added 30 problems
% 041117 med  xxx_prob removed and code added
% 070221 hkh  Change IntVars to use full vector format
% 080603 med  Switched to conAssign, cleaned
% 091015 hkh  Avoid slacks in problem 15 and 16
% 091015 hkh  Avoid expanding KNAPSACK problems to standard form for mipSolve

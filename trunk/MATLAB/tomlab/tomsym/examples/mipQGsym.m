% TomSym version of mipQG

Name='Weingartner 1 - 2/28 0-1 knapsack';

A = [ 45      0     85     150     65     95     30      0    170  0 ...
    40     25     20       0      0     25      0      0     25  0 ...
    165     0     85       0      0      0      0    100  ;
    30     20    125       5     80     25     35     73     12  15 ...
    15     40      5      10     10     12     10      9      0  20 ...
    60     40     50      36     49     40     19    150];
b_U = [600;600];
c   = [1898  440  22507   270  14148  3100  4650  30800   615  4975 ...
    1160 4225    510 11880    479   440   490    330   110   560 ...
    24355 2885  11748  4550    750  3720  1950  10500]'; % 28 weights

toms 28x1 integer x

objective   = -c'*x;
constraints = {A*x<=b_U, 0<=x<=1};

options =  struct;
options.name = Name;
Prob = sym2prob('mip',objective,constraints,[],options);

Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
Prob.Solver.Alg         = 2; % Depth First, then Breadth search
Prob.MIP.KNAPSACK       = 1; % Use knapsack heuristic

Result = tomRun('mipSolve', Prob, 1);
%Result = tomRun('cplex', Prob, 1);
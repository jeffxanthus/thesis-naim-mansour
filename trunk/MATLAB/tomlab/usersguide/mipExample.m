   Name='Weingartner 1 - 2/28 0-1 knapsack';
   % Problem formulated as a minimum problem
   A = [ 45      0     85     150     65     95     30      0    170  0 ...
         40     25     20       0      0     25      0      0     25  0 ...
         165     0     85       0      0      0      0    100  ; ...
         30     20    125       5     80     25     35     73     12  15 ...
         15     40      5      10     10     12     10      9      0  20 ...
         60     40     50      36     49     40     19    150]; 
   b_U = [600;600];  % 2 knapsack capacities
   c   = [1898  440  22507   270  14148  3100  4650  30800   615  4975 ...
          1160 4225    510 11880    479   440   490    330   110   560 ...
         24355 2885  11748  4550    750  3720  1950  10500]'; % 28 weights

   % Make problem on standard form for mipSolve
   [m,n]   = size(A);
   A       = [A eye(m,m)];
   c       = [-c;zeros(m,1)]; % Change sign to make a minimum problem
   [mm nn] = size(A);
   x_L     = zeros(nn,1);
   x_U     = [ones(n,1);b_U];
   x_0     = [zeros(n,1);b_U];

   fprintf('Knapsack problem. Variables %d. Knapsacks %d\n',n,m);
   fprintf('Making standard form with   %d variables\n',nn);

   % All original variables should be integer, but also slacks in this case
   IntVars = nn;  % Could also be set as: IntVars=1:nn; or IntVars=ones(nn,1);
   x_min   = x_L; x_max  = x_U; f_Low  = -1E7; % f_Low <= f_optimal must hold
   n       = length(c);
   b_L     = b_U;
   f_opt   = -141278;


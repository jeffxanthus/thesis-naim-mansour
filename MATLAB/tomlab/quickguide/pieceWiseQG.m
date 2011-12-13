% pieceWiseQG is a small example problem for defining and solving piece-
% wise linear programming problems using TOMLAB. The exmaple was taken from
% the ILOG CPLEX 11.0 User's Manual: a Transport Example.
%
% Objects must be shipped from supply locations to demand locations.
%
% Supply amounts (i): [1000 850 1250]
% Demand amounts (j): [900 1200 600 400]
%
% The model forces all demands to be satisfied and all supplies to be
% shipped.
%
% sum{i=1:nSupply} x(i,j) = supply(i)
% sum{j=1:nDemand} x(i,j) = demand(j)
%
% The cost (to be minimized) is defined as follows:
%
% sum{i=1:nSupply} sum{j=1:nDemand} cost(i,j)*x(i,j)
%
% The cost is a piece-wise linear function decribed by:
%
% Case 1 (non-convex):
% 120 per item, between 0 and 200
% 80 per item, between 200 and 400
% 50 per item, between 400 and above
%
% Case 2 (convex):
% 30 per item, between 0 and 200
% 80 per item, between 200 and 400
% 130 per item, between 400 and above

nDemand = 4;   % Demand nodes
nSupply = 3;   % Supply nodes
supply = [1000;850;1250];
demand = [900;1200;600;400];

N = nDemand*nSupply;
%supply1 -> demand, supply2 -> demand, supply3 -> demand

IntVars = (1:nDemand*nSupply)';

% All supplies shipped to demand centers
A1 = zeros(nSupply,N);
for i=1:nSupply
    A1(i,(i-1)*nDemand+1:i*nDemand) = 1;
end
b_L1 = supply;
b_U1 = supply;

% All supplies received at demand centers
A2 = zeros(nDemand,N);
for i=1:nDemand
    A2(i,i:demandNodes:N) = 1;
end
b_L2 = demand;
b_U2 = demand;

% Merge constraints
A = [A1;A2];
b_L = [b_L1;b_L2];
b_U = [b_U1;b_U2];

% Bounds on variables based on supply
x_L = zeros(N,1);
x_U = repmat(demand,nSupply,1);

% Also limited by demand
for i=1:nSupply
    idx = (i-1)*nDemand+1:i*nDemand;
    x_U(idx) = min(supply(i),x_U(idx));
end

% We need N extra variables. One for each piece-wise function.
% The sum of these is the total cost.
c = [zeros(N,1);ones(N,1)];
A = [A, zeros(size(A,1),N)];
x_L = [x_L; zeros(N,1)];
x_U = [x_U; inf*ones(N,1)];
x_0 = []; %No starting point possible
IntVars = [IntVars; zeros(N,1)];

Prob1 = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'PW-EX1', [], [], IntVars);
Prob2 = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'PW-EX2', [], [], IntVars);
Prob3 = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'PW-EX3', [], [], IntVars);
Prob4 = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'PW-EX4', [], [], IntVars);

% Syntax 1, non-convex
type = 'cplex';
point = [200;400];
slope = [120;80;50];
for i=1:N
    var    = i;   % Variable in piece-wise linear function
    funVar = N+i; % New variable to act as the new pw function.
    % We only need to specify the value of the function at one point,
    % f(a) = fa. This point and the slopes describe the function completely.
    % fa is zero when no units are shipped.
    a  = 0;
    fa = 0;
    Prob1 = addPwLinFun(Prob1, 1, type, var, funVar, point, slope, a, fa);
end

% Solve problem 1 (non-convex, syntax 1)
Result1 = tomRun('cplex', Prob1, 1);

% Syntax 2, non-convex
type = 'cplex';
firstSlope = 120;
point = [200;400];
value = [200*120;200*120+200*80];
lastSlope = 50;
for i=1:N
    var    = i;   % Variable in piece-wise linear function
    funVar = N+i; % New variable to act as the new pw function.
    Prob2 = addPwLinFun(Prob2, 2, type, var, funVar, firstSlope, point, ...
        value, lastSlope);
end

% Solve problem 1 (non-convex, syntax 2)
Result2 = tomRun('cplex', Prob2, 1);

% Syntax 1, convex
type = 'cplex';
point = [200;400];
slope = [30;80;130];
for i=1:N
    var    = i;
    funVar = N+i;
    a      = 0;
    fa     = 0;
    Prob3 = addPwLinFun(Prob3, 1, type, var, funVar, point, slope, a, fa);
end

% Solve problem 1 (convex, syntax 1)
Result3 = tomRun('cplex', Prob3, 1);

% Syntax 4, convex
type = 'cplex';
firstSlope = 30;
point = [200;400];
value = [200*30;200*30+200*80];
lastSlope = 130;
for i=1:N
    var    = i;
    funVar = N+i;
    a      = 0;
    fa     = 0;
    Prob4 = addPwLinFun(Prob4, 2, type, var, funVar, firstSlope, point, ...
        value, lastSlope);
end

% Solve problem 1 (convex, syntax 2)
Result4 = tomRun('cplex', Prob4, 1);

disp(sprintf('Syntax 1, non-convex results \n'));
disp(reshape(Result1.x_k(1:nDemand*nSupply,1),nDemand,nSupply));
disp(sprintf('Syntax 2, non-convex results \n'));
disp(reshape(Result2.x_k(1:nDemand*nSupply,1),nDemand,nSupply));
disp(sprintf('Syntax 1, convex results \n'));
disp(reshape(Result3.x_k(1:nDemand*nSupply,1),nDemand,nSupply));
disp(sprintf('Syntax 2, convex results \n'));
disp(reshape(Result4.x_k(1:nDemand*nSupply,1),nDemand,nSupply));

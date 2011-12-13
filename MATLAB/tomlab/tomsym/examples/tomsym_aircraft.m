%% Aircraft allocation under uncertain demand.
% TomSym implementation of GAMS Example (AIRCRAF,SEQ=8)
%
% The objective of this model is to allocate aircrafts to routes to
% maximize the expected profit when traffic demand is uncertain. Two
% different formulations are used, the delta and the lambda formulation.
%
% Dantzig, G B, Chapter 28. In Linear Programming and Extensions.
% Princeton University Press, Princeton, New Jersey, 1963.
%
% i: aircraft types and unassigned passengers (a, b, c, d)
%
% j: assigned and unassigned routes (1,2,3,4,5)
%
% h: demand states (1,2,3,4,5)

% Parameters

% Demand distribution. Row correspond to route (j), columns to demand
% states (h)
dd = [200     220    250    270    300;
    50     150      0      0      0;
    140     160    180    200    220;
    10      50     80    100    340;
    580     600    620      0      0];

% Probability of demand state (h) for each route (j)
lambda = [.2     .05    .35    .2     .2;
    .3     .7      0      0      0;
    .1     .2     .4     .2     .1;
    .2     .2     .3     .2     .1;
    .1     .8     .1      0      0];

% Aircraft costs (in thousands). The rows are aircraft types and unassigned
% passengers (i). The columns are assigned and unassigned routes (j).
c = [18         21         18          16           10;
    0         15         16          14            9;
    0         10          0           9            6;
    17         16         17          15           10];

% Passenger capacity of aircraft i on route j
p = [16         15          28          23          81;
    0         10          14          15          57;
    0          5           0           7          29;
    9         11          22          17          55];

% Aircraft availability on i
aa = [10;19;25;15];

% Revenue lost (1000 per 100  bumped)
k = [13;13;7;7;1];

% Expected demand
ed = sum(lambda.*dd,2);

% Probability of exceeding demand increment h on route j
g = zeros(size(lambda));
for j = 1:5
    g(:,j) = sum(lambda(:,j:end),2);
end

% Incremental passenger load in demand states
deltb = zeros(size(dd));
for i=1:5
    for j=1:5
        if j > 1
            if dd(i,j)-dd(i,j-1) >= 0
                deltb(i,j) = dd(i,j)-dd(i,j-1);
            end
        else
            deltb(i,j) = dd(i,j);
        end
    end
end

% Variables (positive)

% Number of aircraft type i assigned to route j
toms 4x5 x

% y = passengers actually carried
% b = passengers bumped
toms 5x5 y b

% All variables have to be non-zero.
cbnd = {x >= 0; y >= 0; b >= 0};

% Certain variables have to be zero
cbnd = {cbnd; b(2,3:5) <= 0; b(5,4:5) <= 0}; %From matrix dd
cbnd = {cbnd; y(2,3:5) <= 0; y(5,4:5) <= 0}; %From matrix dd
cbnd = {cbnd; x(2:3,1) <= 0; x(3,3) <= 0};   %From matrix c

% Aircraft balance
eq1 = {sum(x,2) <= aa};

% Demand balance
eq2 = {sum(p.*x,1) >= sum(y,2)'};

% Definition of boarded passengers
eq3 = {y <= repmat(sum(p.*x)',1,5)};

% Definition of bumped passengers
eq4 = {b == dd - y};

% Operating cost definition
oc = sum(sum(c.*x));

% Bumping cost definition: version 1
bc1 = sum(k.*ed) - sum(k.*sum(g.*y,2));

% Bumping cost definition: version 2
bc2 = sum(sum(repmat(k,1,5).*lambda.*b));

% Objective functions (versions 1 and 2)
obj1 = oc+bc1;
obj2 = oc+bc2;

% Version 1 requires upper bounds on y
cbnd1 = {y <= deltb};

% Model 1
solution1 = ezsolve(obj1,{cbnd; cbnd1; eq1; eq2});

% Model 2
solution2 = ezsolve(obj2,{cbnd; eq1; eq3; eq4});
%% A Transportation Problem
% TomSym implementation of GAMS Example (TRNSPORT,SEQ=1)
%
% This problem finds a least cost shipping schedule that meets
% requirements at markets and supplies at factories.
%
% Dantzig, G B, Chapter 3.3. In Linear Programming and Extensions.
% Princeton University Press, Princeton, New Jersey, 1963.
%
% This formulation is described in detail in:
% Rosenthal, R E, Chapter 2: A GAMS Tutorial. In GAMS: A User's Guide.
% The Scientific Press, Redwood City, California, 1988.

% Capacity of plant i in cases
a = [350;600];

% Demand at market j in cases
b = [325;300;275];

% Distance in thousands of miles
d = [2.5 1.7 1.8
    2.5 1.8 1.4];

% Freight in dollars per case per thousand miles
f = 90;

% Transport cost in thousands of dollars per case
c = f*d/1000;

% Shipment quantities in cases
toms 2x3 x
cbnd = (x >= 0);

% Define objective function
cost = sum(sum(c.*x));

% Observe supply limit at plant i
eq1 = {sum(x,2) <= a};

% Satisfy demand at market j
eq2 = {sum(x,1)' >= b};

solution = ezsolve(cost,{cbnd, eq1, eq2});

disp(' ');
disp('Shipment quantities: ');
disp(solution.x);
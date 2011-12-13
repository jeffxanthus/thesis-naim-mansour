%% Area of Hexagon Test Problem
% TomSym implementation of GAMS Example (HIMMEL16,SEQ=36)
%
% The physical problem is to maximize the area of a hexagon in which the
% diameter must be less than or equal to one. The formulation in Himmelblau
% is different from the one given here, because certain fixed variables
% have been eliminated. However, the formulation given here is more natural
% and easier to understand.  It makes use of the fact that one vertex is
% fixed at the origin and uses the vector scalar product to calculate the
% areas of all triangles originating from the origin. All terms in x(1) and
% y(1) thus vanish when the algebraic expression is simplified.
%
% The problem appears many other places, e.g. as example 108 in W. Hock and
% K. Schittkowski: Test Examples for Nonlinear Programming Codes, Lecture
% Notes in Economics and Mathematical Systems, 187, Springer Verlag, 1981,
% and as the test example in P. E. Gill, W. Murray, M. A. Saunders, and
% M. Wright: User's Guide for SOL/NPSOL: A FORTRAN Package for Nonlinear
% Programming, Tech. Rep. 83-12, Dept. of Operation Research, Stanford
% University.
%
% Himmelblau, D M, Problem Number 16. In Applied Nonlinear Programming.
% Mc Graw Hill, New York, 1972.
%
% i: indices for the 6 points (1-6);

% x-coordinates of the points
x    = tom('x',6,1);

% y-coordinates of the points
y    = tom('y',6,1);

% area of the i'th triangle (circular)
trianglearea = tom('trianglearea',6,1);

% Maximal distance between i and j
eq1 = {};
for i=1:6
    for j=1:6
        if i<j
            eq1 = {eq1; (x(i)-x(j))^2+(y(i)-y(j))^2 <= 1};
        end
    end
end

% Area definition for triangle i
eq2 = {trianglearea == 0.5*(x.*y([2:end,1],1)-y.*x([2:end,1],1))};

% Total area
totarea = sum(trianglearea);

% Initial conditions
cbnd = {x(1) == 0; y(1) == 0; y(2) == 0; trianglearea >= 0};

% Starting point
x0 = {x == [0;0.5;0.5;0.5;0;0]
    y == [0;0;0.4;0.8;0.8;0.4]
    trianglearea == 0.5*(x.*y([2:end,1],1)-y(i)*x([2:end,1],1))};

solution = ezsolve(-totarea,{cbnd, eq1, eq2},x0);
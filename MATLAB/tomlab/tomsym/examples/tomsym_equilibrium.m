%% Chemical Equilibrium Problem
% TomSym implementation of GAMS Example (WALL,SEQ=76)
%
% A Sample Nonlinear system to solve Chemical Equilibrium models.
%
% Wall, T W, Greening, D, and Woolsey, R E D, Solving Complex Chemical
% Equilibria Using a Geometric-Programming Based Technique. OR 34, 3
% (1987).

toms ba so4 baoh oh hso4 h

r1 = {ba*so4 == 1};
r2 = {baoh/ba/oh == 4.8};
r3 = {hso4/so4/h == .98};
r4 = {h*oh == 1};

b1 = {ba + 1e-7*baoh == so4 + 1e-5*hso4};
b2 = {2*ba + 1e-7*baoh + 1e-2*h == 2*so4 + 1e-5*hso4 + 1e-2*oh};

x0 = {1 == ba; 1 == so4; 1 == baoh
    1 == oh; 1 == hso4; 1 == h };

solution = ezsolve(ba,{r1,r2,r3,r4,b1,b2},x0);
%% Ramsey Model of Optimal Economic Growth
% TomSym implementation of GAMS Example (RAMSEY,SEQ=63)
%
% This formulation is described in GAMS/MINOS: Three examples by Alan S.
% Manne, Department of Operations Research, Stanford University, May 1986.
%
% Ramsey, F P, A Mathematical Theory of Saving. Economics Journal (1928).
%
% Murtagh, B, and Saunders, M A, A Projected Lagrangian Algorithm and
% its Implementation for Sparse Nonlinear Constraints. Mathematical
% Programming Study 16 (1982), 84-117.
%
% The optimal objective value is 2.4875
%
% t: time periods

time   = (1:11)';
tfirst = time(1);
tlast  = time(end);

bet = 0.95; % discount factor
b = 0.25;   % capital's value share
g = 0.03;   % labor growth rate
ac = 0.15;  % absorptive capacity rate
k0 = 3;     % initial capital
i0 = 0.05;  % initial investment
c0 = 0.95;  % initial consumption

% Discount factor
betavec = bet.^(1:11)';
betavec(end) = betavec(end)/(1-bet);

% The last period's utility discount factor, is calculated by summing the
% infinite geometric series from the horizon date onward. Because of the
% logarithmic form of the utility function, the post-horizon consumption
% growth term may be dropped from the maximand.

% Output-labor scaling vector
a     = (c0+i0)/k0^b;
alt   = a*(1+g).^((1-b)*((1:11)'-1));

% Capital stock, consumption and investment
toms 11x1 kt ct it

% Capacity constraint
eq1 = {alt.*kt.^b  ==  ct + it};

% Capital balance
eq2 = {kt(2:end) == kt(1:end-1)+it(1:end-1)};

% Terminal condition (provides for post-terminal growth)
eq3 = {g*kt(end) <= it(end)};

% Objective
obj = sum(betavec.*log(ct));

% Bounds
cbnd = {kt >= k0; ct >= c0; it >= i0
    it <= i0*((1+ac).^(time-1))};

cbndinit = {kt(1) == k0};

solution = ezsolve(-obj,{eq1, eq2, eq3, cbnd, cbndinit});
%% Optimal Pricing and Extraction for OPEC
% TomSym implementation of GAMS Example (PINDYCK,SEQ=28)
%
% This model finds the optimal pricing and extraction of oil for the OPEC
% cartel.
%
% Pindyck, R S, Gains to Producers from the Cartelization of
% Exhaustible Resources. Review of Economics and Statistics 60
% (1978), 238-251.

t  = (1:17)';

% demand(t) equilibrium world demand for fixed prices
demand = 1+2.3*1.015.^(t-1);

% Due to a bug in the Matlab syntax, the parser cannot know if f(x) is a
% function call or the x:th element of the vector f. So it has to guess.
% The Matlab parser doesn't understand that "toms" creates variables, so it
% may get confused if one of the names is previously used by a function or
% script.  By assigning something to each variable before calling toms, we
% inform the Matlab parser that it is indeed dealing with variables.

p = []; td = []; s = []; cums = []; d = []; r = []; rev = [];

% p(t):   world price of oil
% td(t):  total demand for oil
% s(t):   supply of oil by non-opec countries
% cums(t):  cumulative supply by non-opec countries
% d(t):   demand for opec-oil
% r(t):   opec reserves
% rev(t): revenues in each period
toms 17x1 p td s cums d r rev

% Positive variables
cbnd = {p >= 0; td >= 0; s >= 0; cums >= 0
    d >= 0; r >= 0};

i = 2:17;

% Total demand equation
eq1 =  td(i) == 0.87*td(i-1) - 0.13*p(i) + demand(i-1);

% Supply equation for non-opec countries
eq2 =   s(i) == 0.75*s(i-1) + (1.1+0.1*p(i)).*1.02.^(-cums(i)./7);

% Accounting equation for cumulative supply
eq3 =  cums(i) == cums(i-1) + s(i);

% Accounting equation for opec reserves
eq4 =   r(i) == r(i-1) - d(i);

% Demand equation for opec
eq5 = d(i) == td(i) - s(i);

% Yearly objective function value
eq6 = rev(i) == d(i).*(p(i)-250./r(i));

profit = sum(rev(2:17).*1.05.^(1-(1:16)'));

% Fixed initial conditions
cbnd1 = {td(1) == 18; s(1) == 6.5; p(1) == 0
    r(1) == 500; cums(1) == 0; d(1) == 11.5; rev(1) == 0};

% Starting point for optimization
x0 = {td == 18; s == [6.5;7*ones(16,1)]; cums == 7*(0:16)'
    d == td-s; p == [0;14*ones(16,1)]};

r0 = [500;zeros(16,1)];
for i=2:17
    r0(i) = r0(i-1)-(18-7);
end

i = 2:17;
x0 = {x0; r == r0; rev == [0;d(i).*(p(i)-250./r(i))] };

cons = {cbnd; cbnd1; eq1; eq2; eq3; eq4; eq5; eq6};
solution = ezsolve(-profit,cons,x0);
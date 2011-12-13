%% Maximum Likelihood Estimation
% TomSym implementation of GAMS Example (LIKE,SEQ=25)
%
% This application from the biomedical area tests the hypothesis
% that a population of systolic blood pressure can be separated into
% three distinct groups.
%
% Bracken, J, and McCormick, G P, Chapter 8.5. In Selected Applications of
% Nonlinear Programming. John Wiley and Sons, New York, 1968, pp. 90-92.
%
% i: Observations (1-31)
%
% g: Groups (one, two, three)

% Systolic blood pressure data
pressure = [95 105 110 115 120 125 130 135 140 145 150 155 160 165 170 ...
    175 180 185 190 195 200 205 210 215 220 225 230 235 240 245 260]';
frequency = [1   1   4   4  15  15  15  13  21  12  17   4  20   8  17 ...
    8   6   6   7   4   3   3   8   1   6   0   5   1   7   1   2]';

y = pressure;
w = frequency;

% Constant
c = 1/sqrt(2*3.14159);

% p(g): proportion of population
% m(g): population mean
% s(g): population standard deviation

toms 3x1 p m s

% Maximum likelihood function
toms i
mlf = fsum(lookup(w,i)*log(c*sum(p./s.*exp(-.5*((lookup(y,i)-m)./s).^2))),...
    i, 1:31);

eq1 = {sum(p) == 1};

eq2 = {m(2) >= m(1); m(3) >= m(2)};

eq3 = {p >= 0.1; s >= 0.1; m >= 0};

x0 = {p == 1/3; m == 100+30*(1:3)'; s == 15};

options = struct;
options.solver = 'conopt';
solution = ezsolve(-mlf,{eq1,eq2,eq3},x0,options);
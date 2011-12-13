%% Structural Optimization of Process Flowsheets
% TomSym implementation of GAMS Example (PROCSEL,SEQ=116)
%
% The goal is the profitable production of chemical C, which can be
% produced from chemical B where B may be the raw material that can be
% purchased from the external market or an intermediate that is produced
% from raw material A. There are two alternative paths of producing B from
% A. A mixed-integer nonlinear formulation is presented to solve the
% optimal production and capacity expansion problem.
%
% Kocis, G R, and Grossmann, I E, Relaxation Strategy for the Structural
% Optimization of Process Flow Sheets. Independent Engineering Chemical
% Research 26, 9 (1987), 1869-1880.
%
% Morari, M, and Grossmann, I E, Eds, Chemical Engineering Optimization
% Models with GAMS. Computer Aids for Chemical Engineering Corporation,
% 1991.
%
% Process flowsheet:
%
%          A2    +-----+  B2      BP
%         +----->|  2  |----->+    |
%    A    |      +-----+      |    |  B1    +-----+    C1
%    ---->|                   +----+------->|  1  |-------->
%         |      +-----+      |             +-----+
%         +----->|  3  |----->+
%          A3    +-----+  B3
%
% a2: consumption of chemical a in process 2
% a3: consumption of chemical a in process 3
% b2: production capacity of chemical b in process 2
% b3: production capacity of chemical b in process 3
% bp: amount of chemical b purchased in external market
% b1: consumption of chemical b in process 1
% c1: production capacity of chemical c in process 1

toms a2 a3 b2 b3 bp b1 c1
posvbls = [a2 a3 b2 b3 bp b1 c1]';

% Variables are positive
cbnd1 = {posvbls >= 0};

% Binary variables, existence of process 1 2 and 3
toms integer y1 y2 y3

eq1 = {c1 == 0.9*b1        % Input-output for process 1
    exp(b2) - 1 ==  a2     % Input-output for process 2
    exp(b3/1.2) - 1 == a3  % Input-output for process 3
    b1 == b2 + b3 + bp     % Mass balance for chemical b
    c1  <= 2*y1            % Logical constraint for process 1
    b2  <= 4*y2            % Logical constraint for process 2
    b3  <= 5*y3};          % Logical constraint for process 3

% total profit in million dollars per year
% sales revenue, fixed investment cost, operating cost and purchases
pr = 11*c1;
pr = pr - 3.5*y1 - y2 - 1.5*y3;
pr = pr - b2 - 1.2*b3;
pr = pr - 1.8*(a2+a3) - 7*bp;

% c1 is less than 1
cbnd2 = {c1 <= 1};

solution = ezsolve(-pr,{cbnd1,cbnd2,eq1});
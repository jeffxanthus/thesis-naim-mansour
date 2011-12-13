%% Ajax Paper Company Production Schedule
% TomSym implementation of GAMS Example (AJAX,SEQ=60)
%
% This sample model is taken from the cybernet pds/apex sample library of
% models. A paper manufacturer can produce four different types of paper on
% three different machines. Given a fixed demand schedule the objective is
% to find a production plan that maximizes monthly profits.
%
% CDC, PDS/APEX Sample Model Library, 1977. Control Data Corporation

% m: machines at mills (machine 1, machine 2, machine 3)
%
% g: paper grades (20-bond-wt, 25-bond-wt, c-bond-ext, tissue-wrp)

toms 3x1 m
toms 4x1 g

% Matrix prate(g,m) production rate (tons per hour)
prate = [53        52         49;
    51        49         44;
    52        45         47;
    42        44         40];

% Matrix pcost(g,m) production cost ($ per ton)
pcost = [76        75        73;
    82        80        78;
    96        95        92;
    72        71        70];

demand = [30000; 20000; 12000; 8000];
price = [77;81;99;105];

% Available machine time
avail = [672; 600; 480];

% Production (tons per month)
toms 4x3 outp

% Variable are positive
cbnd = {outp >= 0};

% Machine capacity (hours per month)
eq1 = {sum(outp./prate,1) <= avail'};

% Demand (tons per month)
eq2 = {sum(outp,2) == demand};

% Profit definition ($ per month)
profit = sum(sum(demand.*price)) - sum(sum(pcost.*outp));

solution = ezsolve(-profit,{cbnd, eq1, eq2});
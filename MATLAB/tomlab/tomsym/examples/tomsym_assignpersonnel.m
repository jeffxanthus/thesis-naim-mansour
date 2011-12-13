%% Assigning Personnel to Machines
%
%% Problem description
% An operator needs to be assigned to each of the six machines in a
% workshop. Six workers have been pre-selected. Everyone has
% undergone a test of her productivity on every machine. The table
% below lists the productivities in pieces per hour. The machines run
% in parallel, that is, the total productivity of the workshop is the
% sum of the productivities of the people assigned to the machines.
%
% Productivity in pieces per hour
%
%  +-------+-----------------+
%  |       |    Machines     |
%  +-------+--+--+--+--+--+--+
%  |Workers| 1| 2| 3| 4| 5| 6|
%  +-------+--+--+--+--+--+--+
%  |   1   |13|24|31|19|40|29|
%  |   2   |18|25|30|15|43|22|
%  |   3   |20|20|27|25|34|33|
%  |   4   |23|26|28|18|37|30|
%  |   5   |28|33|34|17|38|20|
%  |   6   |19|36|25|27|45|24|
%  +-------+--+--+--+--+--+--+
%
% The objective is to determine an assignment of workers to machines
% that maximizes the total productivity. We may start by calculating
% a (non-optimal) heuristic solution using the following fairly
% natural method: choose the assignment p -> m with the highest
% productivity, cross out the line p and the column m (since the
% person has been placed and the machine has an operator), and
% restart this process until we have assigned all persons. The
% problem should then be solved to optimality using Mathematical
% Programming. And finally, solve the same problem to optimality,
% but for machines working in series.
%
%% Variables
%
%  prodmat                    The productivity matrix
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Oct 7, 2005.   Last modified Mar 8, 2009.

%% Problem setup
prodmat = [ 13 24 31 19 40 29;...
    18 25 30 15 43 22;...
    20 20 27 25 34 33;...
    23 26 28 18 37 30;...
    28 33 34 17 38 20;...
    19 36 25 27 45 24];

p = size(prodmat,1);  %workers
m = size(prodmat,2);  %machines

assign = tom('assign',p,m,'int');

% All variables are binary
bnds = {0 <= assign <= 1};

% Worker/machine constraints
con = {sum(assign,1) == 1, sum(assign,2) == 1};

% Objective
objective = -sum(sum(prodmat.*assign));

constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Assigning Personnel to Machines';
sol1 = ezsolve(objective,constraints,[],options);
f_k1 = subs(objective,sol1);

% Series
pmin = tom('pmin',1,1);
bnds = {bnds, pmin >= 0};

% Productivity bounds 1
con2 = {sum(prodmat.*assign,2) >= pmin};
constraints = {constraints, con2};
objective = -pmin;
sol2 = ezsolve(objective,constraints,[],options);
f_k2 = subs(objective,sol2);

PriLev = 1;
if PriLev > 0
    m   = size(prodmat,1); % number of machines
    w   = m;               % number of workers
    x1  = sol1.assign;
    disp(['Best parallel work (' num2str(-f_k1) ') when '])
    [workers,machine] = find(x1);
    for i = 1:length(workers)
        disp(['   worker '           num2str(workers(i)) ...
            ' operates machine ' num2str(machine(i))])
    end
    x2  = sol2.assign;
    disp(['Best serial work (' num2str(-f_k2) ') when '])
    [workers,machine] = find(x2);
    for i = 1:length(workers)
        disp(['   worker '           num2str(workers(i)) ...
            ' operates machine ' num2str(machine(i))])
    end
end

% MODIFICATION LOG
%
% 051202 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
%% Planning the Production of Electronic Components
%
%% Problem description
% To augment its competitiveness a small business wishes to improve
% the production of its best selling products. One of its main
% activities is the production of cards with microchips and
% electronic badges. The company also produces the components for
% these cards and badges. Good planning of the production of these
% components therefore constitutes a decisive success factor for the
% company. The demand for components is internal in this case and
% hence easy to anticipate.
%
% For the next six months the production of four products with
% references X43-M1, X43-M2, Y54-N1, Y54-N2 is to be planned. The
% production of these components is sensitive to variations of the
% level of production, and every change leads to a non-negligible
% cost through controls and readjustments. The company therefore
% wishes to minimize the cost associated with these changes whilst
% also taking into account the production and storage costs.
%
% The demand data per time period, the production and storage costs,
% and the initial and desired final stock levels for every product
% are listed in the following table. When the production level
% changes, readjustments of the machines and controls have to be
% carried out for the current month. The cost incurred is
% proportional to the increase or reduction of the production
% compared to the preceding month. The cost for an increase of the
% production is $ 1 per unit but only $ 0.50 for a decrease of the
% production level.
%
% Data for the four products
%
%  +------------------------------------+------------------+-------------+
%  |           Demands                  |        Cost      |     Stock   |
%  +------+----+----+----+----+----+----+----------+-------+-------+-----+
%  | Month|  1 |  2 |  3 |  4 |  5 |  6 |Production|Storage|Initial|Final|
%  +------+----+----+----+----+----+----+----------+-------+-------+-----+
%  |X43-M1|1500|3000|2000|4000|2000|2500|    20    |  0.4  |  10   | 50  |
%  |X43-M2|1300| 800| 800|1000|1100| 900|    25    |  0.5  |   0   | 10  |
%  |Y54-N1|2200|1500|2900|1800|1200|2100|    10    |  0.3  |   0   | 10  |
%  |Y54-N2|1400|1600|1500|1000|1100|1200|    15    |  0.3  |   0   | 10  |
%  +------+----+----+----+----+----+----+----------+-------+-------+-----+
%
% What is the production plan that minimizes the sum of costs
% incurred through changes of the production level, production
% and storage costs?
%
%% Variables
%
%  demand                      the demand for each component and month
%  prodcost                    Cost to produce a component
%  storagecost                 Cost to store a component
%  initialstock                Initial stock
%  finalstock                  Final stock
%  increasecost, decreasecost  Increase or decrease in cost
%                              when producing more or less this month
%                              than last month.
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
demand = [1500 3000 2000 4000 2000 2500;...
    1300  800  800 1000 1100  900;...
    2200 1500 2900 1800 1200 2100;...
    1400 1600 1500 1000 1100 1200];

prodcost        = [20;25;10;15];
storagecost     = [.4;.5;.3;.3];
initialstock    = [10;0;50;0];
finalstock      = [50;10;30;10];
increasecost    = 1;
decreasecost    = 0.5;

t = size(demand,2);
p = size(demand,1);

produce = tom('produce',p,t,'int');
store   = tom('store',p,t,'int');
add     = tom('add',t,1,'int');
reduce  = tom('reduce',t,1,'int');

% All slots are integers
bnds = {produce >= 0, store >= 0, add >= 0, reduce >= 0};
bnds1 = {store(:,end) >= finalstock};

% Production equilibrium constraint at start
con1 = {store(:,1) == initialstock + produce(:,1) - demand(:,1)};

% Production equilibrium in process
con2 = {store(:,2:end) == store(:,1:end-1) + produce(:,2:end) - demand(:,2:end)};

% Add/reduction constraint
con3 = {sum(produce(:,2:end),1)' - sum(produce(:,1:end-1),1)' == ...
    add(2:end) - reduce(2:end)};

% Objective
objective = sum(prodcost'*produce + storagecost'*store) + ...
    sum(increasecost*add(2:end)' + decreasecost*reduce(2:end)');

constraints = {bnds, bnds1, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Production of Electronic Components';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    for m = 1:t,
        disp(['Solution for month ' num2str(m)])
        disp(['produce  ' num2str(sol.produce(:,m)')])
        disp(['stock    ' num2str(sol.store(:,m)')])
        disp(['increase ' num2str(sol.add(m))])
        disp(['decrease ' num2str(sol.reduce(m))])
        disp(' ')
    end
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
% 060203 med   Removed printing of temp
% 090407 med   Converted to tomSym
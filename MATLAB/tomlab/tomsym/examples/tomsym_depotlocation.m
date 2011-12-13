%% Depot Location
%
%% Problem description
% A large company wishes to open new depots to deliver to its sales
% centers. Every new set-up of a depot has a fixed cost. Goods are
% delivered from a depot to the sales centers close to the site.
% Every delivery has a cost that depends on the distance covered. The
% two sorts of cost are quite different: set-up costs are capital
% costs which may usually be written off over several years, and
% transport costs are operating costs. A detailed discussion of how
% to combine these two costs is beyond the scope of this text — we
% assume here that they have been put on some comparable basis,
% perhaps by taking the costs over a year.
%
% There are 12 sites available for the construction of new depots and
% 12 sales centers need to receive deliveries from these depots.
%
% The following table gives the costs (in thousand $) of satisfying
% the entire demand of each customer (sales center) from a depot (not
% the unit costs). So, for instance, the cost per unit of supplying
% customer 9 (who has a total demand of 30 tonnes according to the
% next table) from depot 1 is $ 60000/30t, i.e. $ 2000/t. Certain
% deliveries that are impossible are marked with "Inf".
%
% Delivery costs for satisfying entire demand of customers
%
%  +-----+-----------------------------------------------+
%  |     |               Customer                        |
%  |Depot|  1   2   3   4   5   6   7   8   9  10  11  12|
%  +-----+---+---+---+---+---+---+---+---+---+---+---+---+
%  |   1 |100| 80| 50| 50| 60|100|120| 90| 60| 70| 65|110|
%  |   2 |120| 90| 60| 70| 65|110|140|110| 80| 80| 75|130|
%  |   3 |140|110| 80| 80| 75|130|160|125|100|100| 80|150|
%  |   4 |160|125|100|100| 80|150|190|150|130|Inf|Inf|Inf|
%  |   5 |190|150|130|Inf|Inf|Inf|200|180|150|Inf|Inf|Inf|
%  |   6 |200|180|150|Inf|Inf|Inf|100| 80| 50| 50| 60|100|
%  |   7 |100| 80| 50| 50| 60|100|120| 90| 60| 70| 65|110|
%  |   8 |120| 90| 60| 70| 65|110|140|110| 80| 80| 75|130|
%  |   9 |140|110| 80| 80| 75|130|160|125|100|100| 80|150|
%  |  10 |160|125|100|100| 80|150|190|150|130|Inf|Inf|Inf|
%  |  11 |190|150|130|Inf|Inf|Inf|200|180|150|Inf|Inf|Inf|
%  |  12 |200|180|150|Inf|Inf|Inf|100| 80| 50| 50| 60|100|
%  +-----+---+---+---+---+---+---+---+---+---+---+---+---+
%
% In addition, for every depot we have the following information: the
% fixed cost for constructing the depot that needs to be included
% into the objective function and its capacity limit, all listed in
% the table below.
%
% Fix costs and capacity limits of the depot locations
%
%  +------------+----+----+-----+----+----+----+----+----+----+--
%  |Depot       | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12|
%  +------------+----+----+-----+----+----+----+----+----+----+--
%  |Cost (m$)   |3.5|  9| 10|  4|  3|  9|  9|  3|  4| 10|  9|3.5|
%  |Capacity (t)|300|250|100|180|275|300|200|220|270|250|230|180|
%  +------------+----+----+-----+----+----+----+----+----+----+--
%
% The quantities demanded by the sales centers (customers), are
% summarized in the following table.
%
% Demand data
%
%  +----------+---+--+--+---+---+---+--+--+--+---+--+---+
%  |Customer  |  1| 2| 3|  4|  5|  6| 7| 8| 9| 10|11| 12|
%  +----------+---+--+--+---+---+---+--+--+--+---+--+---+
%  |Demand (t)|120|80|75|100|110|100|90|60|30|150|95|120|
%  +----------+---+--+--+---+---+---+--+--+--+---+--+---+
%
% In every case, the demand of a customer needs to be satisfied but a
% sales center may be delivered to from several depots. Which depots
% should be opened to minimize the total cost of construction and of
% delivery, whilst satisfying all demands?
%
%% Variables
%
%  delcosts                   Costs to deliver to centre from depot
%  idxinf                     Impossible combination of centre/depot
%  delcosts(idxinf)           A large number < Inf
%  buildcosts                 Costs to build a depot=
%  capacity                   Capacities per potential depot
%  demand                     The customers demand in tonnes
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Oct 7, 2005.   Last modified Mar 8, 2009.

%% Problem setup
delcosts     = [ 100  80  50  50  60 100 120  90  60  70  65 110;...
    120  90  60  70  65 110 140 110  80  80  75 130;...
    140 110  80  80  75 130 160 125 100 100  80 150;...
    160 125 100 100  80 150 190 150 130 inf inf inf;...
    190 150 130 inf inf inf 200 180 150 inf inf inf;...
    200 180 150 inf inf inf 100  80  50  50  60 100;...
    100  80  50  50  60 100 120  90  60  70  65 110;...
    120  90  60  70  65 110 140 110  80  80  75 130;...
    140 110  80  80  75 130 160 125 100 100  80 150;...
    160 125 100 100  80 150 190 150 130 inf inf inf;...
    190 150 130 inf inf inf 200 180 150 inf inf inf;...
    200 180 150 inf inf inf 100  80  50  50  60 100];

idxinf = find(delcosts == inf);
delcosts(idxinf) = 1e6;

buildcosts   = [3500 9000 10000 4000 3000 9000 9000 3000 4000 10000 9000 3500]';
capacity     = [ 300  250   100  180  275  300  200  220  270   250  230  180]';
demand       = [ 120   80    75  100  110  100   90   60   30   150   95  120]';

d = length(buildcosts); % DEPOTS
c = length(demand);     % CUSTOMERS

fflow = tom('fflow',d,c);
build = tom('build',d,1,'int');

% Bounds
bnds1 = {0 <= build <= 1};
bnds2 = {0 <= fflow <= 1, fflow(idxinf) <= 0};

% Customer constraint.
con1 = {sum(fflow,1) == 1};

% Capacity constraint.
con2 = {(demand'*fflow')' <= capacity.*build};

% Objective
objective = sum(sum(delcosts.*fflow)) + buildcosts'*build;

constraints = {bnds1, bnds2, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Depot Location';
sol = ezsolve(objective,constraints,[],options);

f_k = subs(objective,sol);
PriLev = 1;
if PriLev > 0
    build   = find(sol.build)'; % the depots to build
    deliver = sol.fflow';
    deliver = deliver(:,build);
    deliver(find(deliver<0.001)) = 0;
    disp(['minimal total cost = ' num2str(f_k) ])
    disp(['build the depots     ' num2str(build)      ])
    [jj,ii] = size(deliver);                 % ii = depots, jj = customers
    for j = 1:jj,
        disp(['let customer ' num2str(j) ' get '])
        for i = 1:ii,
            if deliver(j,i) ~= 0,
                if deliver(j,i) == 1,
                    disp(['   all of its demand from depot ' num2str(build(i)) ])
                else
                    disp(['   ' num2str(deliver(j,i)) ...
                        ' of its demand from depot ' num2str(build(i)) ])
                end
            end
        end
    end
end

% MODIFICATION LOG
%
% 051019 med   Created
% 060112 per   Added documentation
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
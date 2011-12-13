%% Planning the Production of Fiberglass
%
%% Problem description
% A company produces fiberglass by the cubic meter and wishes to plan
% its production for the next six weeks. The production capacity is
% limited, and this limit takes a different value in every time
% period. The weekly demand is known for the entire planning period.
% The production and storage costs also take different values
% depending on the time period. All data are listed in the following
% table.
%
% Data per week
%
%  +----+------------+------+----------+-----------+
%  |    | Production |Demand|Production| Storage   |
%  |Week|capacity(m3)| (m3) |cost($/m3)|cost ($/m3)|
%  +----+------------+------+----------+-----------+
%  |  1 |    140     | 100  |   5      |   0.2     |
%  |  2 |    100     | 120  |   8      |   0.3     |
%  |  3 |    110     | 100  |   6      |   0.2     |
%  |  4 |    100     |  90  |   6      |   0.25    |
%  |  5 |    120     | 120  |   7      |   0.3     |
%  |  6 |    100     | 110  |   6      |   0.4     |
%  +----+------------+------+----------+-----------+
%
% Which is the production plan that minimizes the total cost of
% production and storage?
%
%% Variables
%
%  capacity          Production capacity over time
%  demand            Demand over time
%  prodcost          Cost to produce over time
%  storcost          Cost to store over time
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
capacity = [140;100;110;100;120;100];
demand   = [100;120;100; 90;120;110];
prodcost = [  5;  8;  6;  6;  7;  6];
storcost = [ .2; .3; .2;.25; .3; .4];

% No variables are integers (cost flow system)
m = length(capacity);
flows = tom('flows',m*2-1,1);

% Bounds
bnds = {flows >= 0, flows(1:2:end) <= capacity};

% First node constraint
con1 = {flows(1) - flows(2) == demand(1)};

% Constraints for all other nodes, except final
con2 = {flows(2:2:end-3) + flows(3:2:end-2) - flows(4:2:end-1) == demand(2:end-1)};

% Final node constraint
con3 = {flows(end-1) + flows(end) == demand(end)};

% Objective
objective = prodcost'*flows(1:2:end) + storcost(1:end-1)'*flows(2:2:end-1);

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Production of Fiber Glass';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    Ntimes = length(capacity);
    produce = sol.flows(1:2:end);
    store   = [sol.flows(2:2:end-1); 0];
    for t = 1:Ntimes,
        if produce(t) > 0 | store(t) > 0,
            disp(['during month ' num2str(t) ])
            if produce(t) > 0,
                disp(['   produce  ' num2str(produce(t))])
            end
            if store(t) > 0,
                disp(['   store  ' num2str(store(t)) ' to next month'])
            end
        end
    end
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
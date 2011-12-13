%% Assignment of Production Batches to Machines
%
%% Problem description
% Having determined a set of ten batches to be produced in the next
% period, a production manager is looking for the best assignment of
% these batches to the different machines in his workshop. The five
% available machines in the workshop may each process any of these
% batches, but since they are all models from different years of
% manufacture the speed with which they process the batches changes
% from one machine to another. Furthermore, due to maintenance
% periods and readjustments, every machine may only work a certain
% number of hours in the planning period. The table below lists the
% processing times of the production batches on all machines and the
% capacities of the machines.
%
% Processing durations and capacities (in hours)
%
%  +-------+-----------------------------+--------+
%  |       |    Production batches       |        |
%  |       +--+--+--+--+--+--+--+--+--+--+--------+
%  |Machine| 1| 2| 3| 4| 5| 6| 7| 8| 9|10|Capacity|
%  +-------+--+--+--+--+--+--+--+--+--+--+--------+
%  |   1   | 8|15|14|23| 8|16| 8|25| 9|17|   18   |
%  |   2   |15| 7|23|22|11|11|12|10|17|16|   19   |
%  |   3   |21|20| 6|22|24|10|24| 9|21|14|   25   |
%  |   4   |20|11| 8|14| 9| 5| 6|19|19| 7|   19   |
%  |   5   | 8|13|13|13|10|20|25|16|16|17|   20   |
%  +-------+--+--+--+--+--+--+--+--+--+--+--------+
%
% The production cost of a batch depends on the machine that
% processes it. The hourly usage cost for every machine depends on
% its technology, its age, its consumption of consumables (such as
% electricity, machine oil) and the personnel required to operate it.
% These differences are amplified by the variation in the processing
% durations of a batch depending on the machine. The table below
% lists the production costs in thousand $. On which machine should
% each batch be executed if the production manager wishes to minimize
% the total cost of production?
%
% Production cost depending on the assignment (in k$)
%
%  +-------+-----------------------------+
%  |       |    Production batches       |
%  |       +--+--+--+--+--+--+--+--+--+--+
%  |Machine| 1| 2| 3| 4| 5| 6| 7| 8| 9|10|
%  +-------+--+--+--+--+--+--+--+--+--+--+
%  |   1   |17|21|22|18|24|15|20|18|19|18|
%  |   2   |23|16|21|16|17|16|19|25|18|21|
%  |   3   |16|20|16|25|24|16|17|19|19|18|
%  |   4   |19|19|22|22|20|16|19|17|21|19|
%  |   5   |18|19|15|15|21|25|16|16|23|15|
%  +-------+--+--+--+--+--+--+--+--+--+--+
%
%% Variables
%
%  batches                    Hours needed to produce batches
%  capacity                   Hours available per machine
%  costs                      Cost to produce batches
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Oct 7, 2005.   Last modified Mar 8, 2009.

%% Problem setup
batches         = [ 8  15  14  23   8  16   8  25   9  17;...
    15   7  23  22  11  11  12  10  17  16;...
    21  20   6  22  24  10  24   9  21  14;...
    20  11   8  14   9   5   6  19  19   7;...
    8  13  13  13  10  20  25  16  16  17];

capacity        = [18 ;19 ;25 ;19 ;20];

costs           = [17  21  22  18  24  15  20  18  19  18;...
    23  16  21  16  17  16  19  25  18  21;...
    16  20  16  25  24  16  17  19  19  18;...
    19  19  22  22  20  16  19  17  21  19;...
    18  19  15  15  21  25  16  16  23  15];

m = length(capacity);
p = size(batches,2);

use = tom('use',m,p,'int');
% All variables are binary.
bnds = {0 <= use <= 1};

% Machine constr.
con1 = {sum(use,1) == 1};

% Batch constr.
con2 = {sum(batches.*use,2) <= capacity};

% Objective
objective = sum(sum(costs.*use));

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Production Batches to Machines';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    [m,b] = size(batches);
    temp  = sol.use;
    for i = 1:m,
        tempidx = find(temp(i,:)>0.5);
        disp([' machine ' num2str(i) ' deals with ' num2str(tempidx)])
    end
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060111 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
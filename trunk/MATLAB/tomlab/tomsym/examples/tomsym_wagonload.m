%% Wagon Load Balancing
%
%% Problem description
% Three railway wagons with a carrying capacity of 100 quintals (1
% quintal = 100 kg) have been reserved to transport sixteen boxes.
% The weight of the boxes in quintals is given in the following
% table. How shall the boxes be assigned to the wagons in order to
% keep to the limits on the maximum carrying capacity and to minimize
% the heaviest wagon load?
%
% Weight of boxes
%
%  +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
%  |Box   | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|
%  +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
%  |Weight|34| 6| 8|17|16| 5|13|21|25|31|14|13|33| 9|25|25|
%  +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
%
% Before implementing a Mathematical Programming solution, one may
% wish to try to see whether it is possible to solve this problem
% instance with the following simple heuristic: until all boxes are
% distributed to the wagons we choose in turn the heaviest unassigned
% box and put it onto the wagon with the least load.
%
%% Variables
%
%  minload         Minimal load accepted per wagon
%  maxload         Maximal load accepted per wagon
%  weights         Weight of each box
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
minload = [  0;  0;  0;];
maxload = [100;100;100;];
weights = [34;6;8;17;16;5;13;21;25;31;14;13;33;9;25;25];

w = length(minload);
b = length(weights);

load = tom('load',b,w,'int');
maxweight = tom('maxweight',1,1,'int');

% Bounds
bnds = {0 <= load <= 1, sum(weights)/3 <= maxweight <= 100};

% Load constraint, exactly one box per wagon
con1 = {sum(load,2) == 1};

% Wagon constraint
con2 = {sum(load.*repmat(weights,1,3),1) <= maxweight};

% Objective
objective = maxweight;

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Wagon Load Balancing';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    wagons = length(minload);
    temp   = sol.load;
    for w = 1:wagons,
        idx = find(temp(:,w) > 0.5);
        disp(['wagon ' num2str(w) ...
            ' has a total load of ' num2str(sum(weights(idx))) ...
            ' by carrying the boxes: ' num2str(idx') ])
    end
end

% MODIFICATION LOG
%
% 051007 med   Created.
% 060111 per   Added documentation.
% 060125 per   Moved disp to end
% 060203 med   Printing updated
% 090308 med   Converted to tomSym
% function Result = wagonloadbalancingEx(PriLev)
%
% Creates a TOMLAB LP problem for wagon loading. Minimizes the max weight
%
% WAGON LOAD BALANCING
% 
% Three railway wagons with a carrying capacity of 100 quintals (1
% quintal = 100 kg) have been reserved to transport sixteen boxes.
% The weight of the boxes in quintals is given in the following
% table. How shall the boxes be assigned to the wagons in order to
% keep to the limits on the maximum carrying capacity and to minimize
% the heaviest wagon load?
%
% Weight of boxes
%
% +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
% |Box   | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|16|
% +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
% |Weight|34| 6| 8|17|16| 5|13|21|25|31|14|13|33| 9|25|25|
% +------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
%
% Before implementing a Mathematical Programming solution, one may
% wish to try to see whether it is possible to solve this problem
% instance with the following simple heuristic: until all boxes are
% distributed to the wagons we choose in turn the heaviest unassigned
% box and put it onto the wagon with the least load.
%
%
% VARIABLES
%
% minload                    Minimal load accepted per wagon
% maxload                    Maximal load accepted per wagon
% weights                    Weight of each box
%
% RESULTS
%
% For an interpretation of the results, set PriLev > 1,
% Result = wagonloadbalancingEx(2);
% 
% REFERENCES
%
% Applications of optimization... Gueret, Prins, Seveaux
% http://web.univ-ubs.fr/lester/~sevaux/pl/index.html
%
% INPUT PARAMETERS
% PriLev       Print Level
%
% OUTPUT PARAMETERS
% Result       Result structure.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 7, 2005.   Last modified Feb 3, 2006.

function Result = wagonloadbalancingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

minload    = [  0;  0;  0;];
maxload    = [100;100;100;];
weights    = [34;6;8;17;16;5;13;21;25;31;14;13;33;9;25;25];

Prob = wagonloadbalancing(minload, maxload, weights);
Prob.SolverInf = 'cplex';
Result = infLinSolve(Prob, PriLev);

if PriLev > 1,
   boxes  = length(weights);
   wagons = length(minload);
   temp   = reshape(Result.x_k,boxes,wagons);
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
% function Result = prodplanwithpersonnelassignEx(PriLev)
%
% Creates a TOMLAB MILP problem for production planning with personnel assignment
%
% PRODUCTION PLANNING WITH PERSONNEL ASSIGNMENT
%
% The company Line Production decides to plan the production of four
% of its products (P1, P2, P3, P4) on its five production lines (L1
% to L5). The company gains net profits of $ 7 for the products P1
% and P4, $ 8 for P2, and $ 9 for P3. The maximum times during
% which the five production lines may operate are different during
% the planning period. The maximum capacities for L1 to L5 are 4500
% hours, 5000 hours, 4500 hours, 1500 hours, and 2500 hours 
% respectively. The table below lists the processing time in hours
% necessary for the production of one unit of every product on every
% production line. Which quantities of P1 to P4 should be produced to
% maximize the total profit? If subsequently a transfer of personnel
% (and hence of working hours) is authorized between production lines
% during the planning period as shown in the second table below,
% which is the maximum profit? How many hours are transfered and
% under what conditions?
%
% Unitary processing times
%
% +--------+-------------------+
% |        |   Lines           |
% +--------+---+---+---+---+---+
% |Products| L1| L2| L3| L4| L5|
% +--------+---+---+---+---+---+
% |   P1   |1.3|0.9|2.0|0.3|0.9|
% |   P2   |1.8|1.7|1.4|0.6|1.1|
% |   P3   |1.3|1.2|1.3|1.0|1.4|
% |   P4   |0.9|1.1|1.0|0.9|1.0|
% +--------+---+---+---+---+---+
%
% Possible transfers of personnel
%
% +------+-------------------+------------------+
% |      |  Destination      |                  |
% +------+---+---+---+---+---+Maximum number of |
% |Origin| L1| L2| L3| L4| L5|transferable hours|
% +------+---+---+---+---+---+------------------+
% |  L1  |  -|yes|yes|yes| no|     400          |
% |  L2  | no|  -|yes| no|yes|     800          |
% |  L3  |yes|yes|  -|yes| no|     500          |
% |  L4  | no| no| no|  -|yes|     200          |
% |  L5  |yes|yes|yes| no|  -|     300          |
% +------+---+---+---+---+---+------------------+
% 
% VARIABLES
%
% profit                     Profit per product
% capacity                   Hours at each site
% timemat                    Time to produce products
% transfermat                Transfer matrix
% maxtransfer                Maximum hours to transfer
%
% RESULTS
%
% For an interpretation of the results, run:
% [Result1, Result2] = prodplanwithpersonnelassignEx(2);
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
% Result       Result structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function [Result1, Result2] = prodplanwithpersonnelassignEx(PriLev)

if nargin < 1
   PriLev = 1;
end

profit        = [7 8 9 7]';
capacity      = [4500 5000 4500 1500 2500]';

timemat       = [ 1.3 0.9 2.0 0.3 0.9;...
      1.8 1.7 1.4 0.6 1.1;...
      1.3 1.2 1.3 1.0 1.4;...
      0.9 1.1 1.0 0.9 1.0];

transfermat   = [0 1 1 1 0;...
      0 0 1 0 1;...
      1 1 0 1 0;...
      0 0 0 0 1;...
      1 1 1 0 0];

maxtransfer   = [400 800 200 500 300]';

flag          = 0;

Prob = prodplanwithpersonnelassign(profit, capacity, timemat,...
   transfermat, maxtransfer, flag);
Result1 = tomRun('cplex', Prob, PriLev);

flag          = 1;

Prob = prodplanwithpersonnelassign(profit, capacity, timemat,...
   transfermat, maxtransfer, flag);
Result2 = tomRun('cplex', Prob, PriLev);


if PriLev > 1,
   disp('without transfer')
   x1 = Result1.x_k;
   for i = 1:length(x1),
      if x1(i) > 0,
         disp(['    produce ' num2str(x1(i)) ' units of P' num2str(i)  ])
      end
   end
   disp('with transfer')
   x2 = Result2.x_k;
   length(x2);
   n1    = 5;
   n2    = 4;
   n     = n1 + n2*n2 + n2;   %amount, transfers, hours
   x2 = Result2.x_k;
   for i = 1:length(x1),
      if x2(i) > 0,
         disp(['    produce ' num2str(x2(i)) ' units of P' num2str(i)  ])
      end
   end
   transfer = reshape(x2(n1:n1+n1*n1-1),n1,n1);
   for i = 1:size(transfer,1),
      idx = find(transfer(i,:));
      for j = 1: length(idx),
         disp(['    transfer ' num2str(transfer(i,idx(j))) ...
               ' hours from ' num2str(i) ' to ' num2str(idx(j)) ])
      end
   end
   hours = x2(n1+n1*n1:end);
   for j = 1:length(hours),
      disp(['    work ' num2str(hours(j)) ' hours at L' num2str(j) ])
   end
end

% MODIFICATION LOG
%
% 051205 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end
% 060131 per   Really moved disp to end

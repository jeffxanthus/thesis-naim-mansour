% function Result = bargeloadingEx(PriLev)
%
% Creates a TOMLAB LP problem for barge loading, maximizing the profit
%
% BARGE LOADING
%
% A shipper on the river Rhine owns a barge of carrying capacity 1500
% m3. Over time he has specialized in the transport of wheat. He has
% seven regular customers who load and unload practically at the same
% places. The shipper knows his costs of transport from long
% experience and according to his personal preferences has concluded
% agreements with his clients for the price charged to them for the
% transport of their wheat. The following table summarizes the
% information about the seven clients. Every client wishes to
% transport a certain number of lots, deciding himself the size of
% his lots in m3. The table lists the price charged by the shipper
% for every transported lot. The last column of the table contains
% the cost incurred by the shipper per transported m3. This cost
% differs depending on the distance covered.
%
% Lots to transport
%
% +------+------------------+--------+----------+--------------+
% |      |Available quantity|Lot size|Price per | Transport    |
% |Client| (no. of lots)    |(in m3) |lot (in $)|cost (in $/m3)|
% +------+------------------+--------+----------+--------------+
% |  1   |     12           |   10   |  1000    |     80       |
% |  2   |     31           |    8   |   600    |     70       |
% |  3   |     20           |    6   |   600    |     85       |
% |  4   |     25           |    9   |   800    |     80       |
% |  5   |     50           |   15   |  1200    |     73       |
% |  6   |     40           |   10   |   800    |     70       |
% |  7   |     60           |   12   |  1100    |     80       |
% +------+------------------+--------+----------+--------------+
%
% The objective of the shipper is to maximize his profit from
% transporting the wheat with lots that may be divided.
%
% Question 1:
% As a first step, assuming that every client has an unlimited
% quantity of wheat, which clients’ wheat should be transported?
%
% Question 2:
% If in addition the actual availability of wheat lots from the
% customers is taken into account, which strategy should the shipper
% adopt?
%
% Question 3:
% What happens if the lots cannot be divided?
%
% VARIABLES
%
% capacity                   = 1500;
% units                      = [12;31;20;25;50;40;60];
% sizes                      = [10;8;6;9;15;10;12];
% unitprice                  = [10;6;6;8;12;8;11]*100;
% cost                       = [80;70;85;80;73;70;80];
% var1                       30 plus
% var2                       30 minus
%
% RESULTS
%
% x_k is a long vector of values, and if we run
% Result      = bargeloadingEx(2);
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
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Result = bargeloadingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

capacity    = 1500;
units       = [12;31;20;25;50;40;60];
sizes       = [10;8;6;9;15;10;12];
unitprice   = [10;6;6;8;12;8;11]*100;
cost        = [80;70;85;80;73;70;80];

Prob = bargeloading(capacity, units, sizes, unitprice, cost);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   disp('Best buy')
   for i = 1:length(units),
      if Result.x_k(i) > 0,
         disp([' ' num2str(Result.x_k(i)) ' unit(s) from ' num2str(i)])
      end
   end
end

% MODIFICATION LOG
%
% 051007 med   Created.
% 060111 per   Added documentation.
% 060125 per   Moved disp to end

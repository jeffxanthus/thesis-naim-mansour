% function Result = airlinehublocationEx(PriLev)
%
% Creates a TOMLAB MIP problem for airline hub location
%
% AIRLINE HUB LOCATION
%
% FAL (the Flying Air Lines) specializes in freight transport. The
% company links the major French cities with cities in the United
% States, namely: Atlanta, Boston, Chicago, Marseille, Nice, and
% Paris. The average quantities in tonnes transported every day by
% this company between these cities are given in the following table.
%
% Average quantity of freight transported between every pair of cities
%
% +---------+-------+------+-------+---------+----+-----+
% |         |Atlanta|Boston|Chicago|Marseille|Nice|Paris|
% |Atlanta  |   0   |  500 | 1000  |   300   |400 |1500 |
% |Boston   |1500   |    0 |  250  |   630   |360 |1140 |
% |Chicago  | 400   |  510 |    0  |   460   |320 | 490 |
% |Marseille| 300   |  600 |  810  |    0    |820 | 310 |
% |Nice     | 400   |  100 |  420  |   730   |  0 | 970 |
% |Paris    | 350   | 1020 |  260  |   580   |380 |   0 |
% +---------+-------+------+-------+---------+----+-----+
%
% We shall assume that the transport cost between two cities i and j
% is proportional to the distance that separates them. The distances
% in miles are given in the next table.
%
% Distances between pairs of cities
%
% +---------+------+-------+---------+----+-----+
% |         |Boston|Chicago|Marseille|Nice|Paris|
% +---------+------+-------+---------+----+-----+
% |Atlanta  |  945 | 605   | 4667    |4749| 4394|
% |Boston   |  866 |3726   | 3806    |3448|     |
% |Chicago  | 4471 |4541   | 4152    |    |     |
% |Marseille|  109 | 415   |         |    |     |
% |Nice     |  431 |       |         |    |     |
% +---------+------+-------+---------+----+-----+
%
% The airline is planning to use two cities as connection platforms
% (hubs) to reduce the transport costs. Every city is then assigned
% to a single hub. The traffic between cities assigned to a given
% hub H1 to the cities assigned to the other hub H2 is all routed
% through the single connection from H1 to H2 which allows the
% airline to reduce the transport cost. We consider that the
% transport cost between the two hubs decreases by 20%. Determine the
% two cities to be chosen as hubs in order to minimize the transport
% cost.
%
% VARIABLES
%
% frights                    Goods between cities
% distance                   Distances
% hubs                       Number of hubs
%
% RESULTS
%
% The result from this example is quite hard to interpret since
% Result.x_k should be split into a tensor of rank 4 with the indices
% i, j, k and l and also into a vector of length c (number of
% cities). The vector indicates the hubs. For all non-zero element
% (i,j,k,l) in the tensor the following interpretation is possible:
% the route from city i to j passes through the hubs k and l.
%
% In the example below temp(3,5,2,6) == 1 since the route from
% 3 (Chicago) to 5 (Nice) pass through 2 (Boston) and 6 (Paris).
% Also temp(4,5,6,6) == 1 since the route from 4 (Marseille) to
% 5 (Nice) pass through only 6 (Paris).
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
% Written Nov 7, 2005.   Last modified Feb 3, 2006.

function Result = airlinehublocationEx(PriLev)

if nargin < 1
   PriLev = 1;
end

frights       = [   0   500   1000    300    400   1500;...
      1500     0    250    630    360   1140;...
      400   510      0    460    320    490;...
      300   600    810      0    820    310;...
      400   100    420    730      0    970;...
      350  1020    260    580    380      0];

distance      = [   0   945   605   4667   4749    4394;...
      945     0   866   3726   3806    3448;...
      605   866     0   4471   4541    4152;...
      4667  3726  4471      0    109     415;...
      4749  3806  4541    109      0     431;...
      4394  3448  4152    415    431       0];

hubs          = 2;

Prob = airlinehublocation(frights, distance, hubs);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1
   c       = 6;   % number of cities
   c_names = ['Atla';'Bost';'Chic';'Mars';'Nice';'Pari'];
   temp    = reshape(Result.x_k(1:c^4),c,c,c,c);
   hubs    = find(Result.x_k(c^4+1:c^4+c));
   
   disp('Put hubs in ')
   for h = 1: length(hubs),
      disp([' ' num2str(c_names(hubs(h),:)) ])
   end
   
   % displaying the routes between cities
   for i = 1:c,
      for j = i+1:c,
         for k = 1:c,
            for l = 1:c,
               if temp(i,j,k,l) == 1,
                  % route involving one hub
                  if k==l & i~=k & j~=k & i~=j,
                     disp(['the route ' num2str([c_names(i,:),' ',   ...
                              c_names(j,:)]) ' requires hub  ' num2str(c_names(k,:))])
                     % route involving two hubs
                  elseif i~=j & i~=k & i~=l & j~=k & j~=l & k~=l,
                     disp(['the route ' num2str([c_names(i,:),' ', ...
                              c_names(j,:)]) ' requires hubs ' ...
                              num2str([c_names(k,:),' ', c_names(l,:)])])
                     % route between hubs
                  elseif ( (i==k & j==l) | (i==l & j==k) ) & k~=l,
                     disp(['the route ' num2str([c_names(i,:),' ', ...
                              c_names(j,:)]) ' is between hubs' ])
                  end
               end
            end
         end
      end
   end
end

% MODIFICATION LOG
%
% 051107 med   Created.
% 060113 per   Added documentation.
% 060116 per   Minor changes.
% 060125 per   Moved disp to end
% 060203 med   Updated print level
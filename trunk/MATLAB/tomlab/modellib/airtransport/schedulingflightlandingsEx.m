% function Result = schedulingflightlandingsEx(PriLev)
%
% Creates a TOMLAB MIP problem for scheduling flight landings
%
% SCHEDULING FLIGHT LANDINGS
%
% The plane movements in large airports are subject to numerous
% security constraints. The problem presented in this section consist
% of calculating a schedule for flight landings on a single runway.
% More general problems have been studied but they are fairly
% complicated (dynamic cases, for instance through planes arriving
% late, instances with several runways, etc.). 
%
% Ten planes are due to arrive. Every plane has an earliest arrival
% time (time when the plane arrives above the zone if traveling at
% maximum speed) and a latest arrival time (influenced among other
% things by its fuel supplies). Within this time window the airlines
% choose a target time, communicated to the public as the flight
% arrival time. The early or late arrival of an aircraft with respect
% to its target time leads to disruption of the airport and causes
% costs. To take into account these cost and to compare them more
% easily, a penalty per minute of early arrival and a second penalty
% per minute of late arrival are associated with every plane. The
% time windows (in minutes from the start of the day) and the
% penalties per plane are given in the following table.
%
% Characteristics of flight time windows
% +-----------------+---+---+---+---+---+---+---+---+---+---+
% |Plane            |  1|  2|  3|  4|  5|  6|  7|  8|  9| 10|
% +-----------------+---+---+---+---+---+---+---+---+---+---+
% |Earliest arrival |129|195| 89| 96|110|120|124|126|135|160|
% |Target time      |155|258| 98|106|123|135|138|140|150|180|
% |Latest Arrival   |559|744|510|521|555|576|577|573|591|657|
% |Earliness penalty| 10| 10| 30| 30| 30| 30| 30| 30| 30| 30|
% |Lateness penalty | 10| 10| 30| 30| 30| 30| 30| 30| 30| 30|
% +-----------------+---+---+---+---+---+---+---+---+---+---+
%
% Due to turbulence and the duration of the time during which a plane
% is on the runway, a security interval has to separate any two
% landings. An entry in line p of column q in the following table
% denotes the minimum time interval (in minutes) that has to lie
% between the landings of planes p and q, even if they are not
% consecutive. Which landing schedule minimizes the total penalty
% subject to arrivals within the given time windows and the required
% intervals separating any two landings?
%
% Matrix of minimum intervals separating landings
%
% +--+--+--+--+--+--+--+--+--+--+--+
% |  | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|
% +--+--+--+--+--+--+--+--+--+--+--+
% | 1| –| 3|15|15|15|15|15|15|15|15|
% | 2| 3| –|15|15|15|15|15|15|15|15|
% | 3|15|15| –| 8| 8| 8| 8| 8| 8| 8|
% | 4|15|15| 8| –| 8| 8| 8| 8| 8| 8|
% | 5|15|15| 8| 8| –| 8| 8| 8| 8| 8|
% | 6|15|15| 8| 8| 8| –| 8| 8| 8| 8|
% | 7|15|15| 8| 8| 8| 8| –| 8| 8| 8|
% | 8|15|15| 8| 8| 8| 8| 8| –| 8| 8|
% | 9|15|15| 8| 8| 8| 8| 8| 8| –| 8|
% |10|15|15| 8| 8| 8| 8| 8| 8| 8| –|
% +--+--+--+--+--+--+--+--+--+--+--+
%
% VARIABLES
%
% costs                      Cost or being late
% mintimes                   Minimal intervals between landings
% var1                       30 plus
% var2                       30 minus
%
% RESULTS
%
% For an interpretation of the results, set PriLev > 1, for example:
% Result = schedulingflightlandingsEx(2);
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
% Written Oct 21, 2005.   Last modified Feb 3, 2006.

function Result = schedulingflightlandingsEx(PriLev)

if nargin < 1
   PriLev = 1;
end

costs       = [ 129 195  89  96 110 120 124 126 135 160;...
      155 258  98 106 123 135 138 140 150 180;...
      559 744 510 521 555 576 577 573 591 657;...
      10  10  30*ones(1,8);...
      10  10  30*ones(1,8)];

mintimes    = [  0  3 15 15 15 15 15 15 15 15;...
      3  0 15 15 15 15 15 15 15 15;...
      15 15  0  8  8  8  8  8  8  8;...
      15 15  8  0  8  8  8  8  8  8;...
      15 15  8  8  0  8  8  8  8  8;...
      15 15  8  8  8  0  8  8  8  8;...
      15 15  8  8  8  8  0  8  8  8;...
      15 15  8  8  8  8  8  0  8  8;...
      15 15  8  8  8  8  8  8  0  8;...
      15 15  8  8  8  8  8  8  8  0];

Prob = schedulingflightlandings(costs, mintimes);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1
   planes       = 10; % number of planes
   [times,idx]  = sort(Result.x_k(1:planes));     
   late         = Result.x_k(planes*planes+2*planes+1:planes*planes+3*planes);
   early        = Result.x_k(planes*planes+1*planes+1:planes*planes+2*planes);
   for p = 1:planes,
      time  = times(p);
      plane = idx(p);
      if (late(plane) > 0.5),
         disp(['plane ' num2str(plane) ' will arrive minute ' ...
               num2str(time) ' (' num2str(late(plane)) ' minute(s) late)'])
      elseif (early(plane) > 0.5),
         disp(['plane ' num2str(plane) ' will arrive minute ' ...
               num2str(time) ' (' num2str(early(plane)) ' minute(s) early)'])
      else
         disp(['plane ' num2str(plane) ' will arrive minute ' ...
               num2str(time) ' (on time)' ])
      end
   end
end

% MODIFICATION LOG
%
% 051021 med   Created.
% 060113 per   Added documentation
% 060125 per   Moved disp to end
% 060203 med   Print level updated
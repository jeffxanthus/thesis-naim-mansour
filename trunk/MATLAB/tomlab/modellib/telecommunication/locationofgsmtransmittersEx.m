% function Result = locationofgsmtransmittersEx(PriLev)
%
% Creates a TOMLAB MIP problem for location of gsm transmitters
%
% LOCATION OF GSM TRANSMITTERS
%
% A mobile phone operator decides to equip a currently uncovered
% geographical zone. The management allocates a budget of $ 10
% million to equip this region. A study shows that only 7 locations
% are possible for the construction of the transmitters and it is
% also known that every transmitter only covers a certain number of
% communities. The figure below represents a schematic map of the
% region with the division into communities and the possible
% locations for transmitters. Every potential site is indicated by 
% four stars and a number, every community is represented by a
% polygon. The number in the center of a polygon is the number of the
% community. 
%
% Certain geographical and topological constraints add to the
% construction cost and reduce the reach of the GSM transmitters. The
% table below lists the communities covered and the cost for every
% site.
% 
% For every community the number of inhabitants is known (see table).
% Where should the transmitters be built to cover the largest
% population with the given budget limit of $ 10M?
%
% Map of the region to cover
%
% +--------+-------+----------+--------+
% |        |       *   7     /         |
% |  1     |   4  *3*       /          |
% |        *       * /--\  /*      11  |
% +-------*1*      |/ 10 \/*6*         |
% |        *\     / \    / -*--+-------+
% |          \   /   \--/      \       |
% |           \ /\     |        \  15  |
% |    2      /   \ 8  |*        \     |
% |          /     \   *5*  12   *\    +
% |         /       \  |*       *7*\  /|
% |        /   5     \ |         *  \/ |
% +-----*-+       *   \|-----+----- /  |
% |    *2* \    -*4*---+     |     \ 14|
% |     *   \  /  *   /      |      \  |
% |          \/      /   9   |       \ |
% |   3       \  6  /        |   13   \|
% |            \   /         |         +
% +-------------+-+----------+---------+
%
% Inhabitants of the communities
%
% +--------------------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
% |Community           | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|
% +--------------------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
% |Population (in 1000)| 2| 4|13| 6| 9| 4| 8|12|10|11| 6|14| 9| 3| 6|
% +--------------------+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
%
% Cost and communities covered for every site
%
% +-------------------+-----+-----+--------+-------+------+-------------+-----------+
% |Site               |  1  |  2  |   3    |   4   |  5   |      6      |     7     |
% +-------------------+-----+-----+--------+-------+------+-------------+-----------+
% |Cost (in million $)| 1.8 | 1.3 |  4.0   |  3.5  | 3.8  |     2.6     |    2.1    |
% |Communities covered|1,2,4|2,3,5|4,7,8,10|5,6,8,9|8,9,12|7,10,11,12,15|12,13,14,15|
% +-------------------+-----+-----+--------+-------+------+-------------+-----------+
%                                                                                  
% VARIABLES
%
% budget                     Money available
% cost                       Cost to build site
% connections                Communities covered by a site
% population                 People living in communities
% var1                       30 plus
% var2                       30 minus
%
% RESULTS
%
% To interpret the results try the following:
% Result   = locationofgsmtransmittersEx(2);
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
% Written Nov 30, 2005.   Last modified Nov 30, 2005.

function Result = locationofgsmtransmittersEx(PriLev)

if nargin < 1
   PriLev = 1;
end

budget      = 10e6;
cost        = [1.8 1.3 4.0 3.5 3.8 2.6 2.1]'*1e6;

connections = [ 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0;...
      0 1 1 0 1 0 0 0 0 0 0 0 0 0 0;...
      0 0 0 1 0 0 1 1 0 1 0 0 0 0 0;...
      0 0 0 0 1 1 0 1 1 0 0 0 0 0 0;...
      0 0 0 0 0 0 0 1 1 0 0 1 0 0 0;...
      0 0 0 0 0 0 1 0 0 1 1 1 0 0 1;...
      0 0 0 0 0 0 0 0 0 0 0 1 1 1 1];

population  = [ 2  4 13  6  9  4  8 12 10 11  6 ...
      14  9  3  6]'*1000;

Prob = locationofgsmtransmitters(budget, cost, connections, ...
   population);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   [loc,com] = size(connections);
   build     = Result.x_k(1:loc);
   cover     = Result.x_k(loc+1:end);
   disp(['Build at locations ' num2str(find(build'))])
   disp(['Cover communities ' num2str(find(cover'))])
end

% MODIFICATION LOG
%
% 051130 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end

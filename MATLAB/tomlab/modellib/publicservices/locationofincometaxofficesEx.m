% function Result = locationofincometaxofficesEx(PriLev)
%
% Creates a TOMLAB MILP problem for location of income tax offices
%
% LOCATION OF INCOME TAX OFFICES
%
% The income tax administration is planning to restructure the
% network of income tax offices in a region. The graph in the figure
% below shows the cities in the region and the major roads. The 
% numbers within () close to the cities indicate the population in
% thousands of inhabitants. The arcs are labeled with the distances
% in kilometers. The income tax administration has determined that
% offices should be established in three cities to provide sufficient
% coverage. Where should these offices be located to minimize the
% average distance per inhabitant to the closest income tax office?
%
% Graph of towns and roads of the region
% 
% 
%   (15)   (10)   (12)   (18)
%   1 -15- 2 -22- 3 -18- 4 
% 
%    | \       / |     |
%    |  24   16  |     |
%   18   \   /   |     |
%    |           20    12
%    |     5(5)  |     |
%    |           |     |
%    |     | \   |     |
%    |    12  24 |     |
%    |     |   \ |     |
%
%   7 -15- 8 -30- 9 -12- 6 (24)
% (11)    (16)   (13)
%    |     |   /  |   /
%    22   25  19  19 22
%    |     | /    | /
%
%  10 -19- 11 -21- 12 (20)
%  (22)    (19)
% 
% VARIABLES
%
% population                 Population of each town
% numloc                     Number of offices to start
% lengths                    The length of the roads
% in/out                     A road i goes between towns
%                            in(i) and out(i)
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = locationofincometaxofficesEx(2);
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
% Written Dec 6, 2005.   Last modified Dec 6, 2005.

function Result = locationofincometaxofficesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

population = [15 10 12 18 5 24 11 16 13 22 19 20]';
numloc     = 3;
in         = [1 1 1 2 3 3 3 4 5 5 6  6 7  7 8  8  9  9 10 11]';
out        = [2 5 7 3 4 5 9 6 8 9 9 12 8 10 9 11 11 12 11 12]';
lengths    = [15 24 18 22 18 16 20 12 12 24 12 22 15 22 30 ...
      25 19 19 19 21]';

Prob = locationofincometaxoffices(population, numloc, in, out, lengths);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   cities = length(population);
   temp   = Result.x_k(1:cities);
   build  = find(temp);
   goto   = reshape(Result.x_k(1+cities:end),cities,cities);
   disp(['Build the offices in towns ' num2str(build') ' and let'])
   for i = 1:length(build),
      disp(['   people from ' num2str(find(goto(build(i),:))) ...
            ' travel to ' num2str(build(i)) ])
   end
end

% MODIFICATION LOG
%
% 051206 med   Created.
% 060118 per   Added documentation.
% 060125 per   Moved disp to end

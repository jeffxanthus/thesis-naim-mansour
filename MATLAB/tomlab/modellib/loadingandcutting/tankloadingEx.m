% function Result = tankloadingEx(PriLev)
%
% Creates a TOMLAB MIP problem for tank loading, maximizing the left over
% capacity
%
% TANK LOADING
%
% Five tanker ships have arrived at a chemical factory. They are
% carrying loads of liquid products that must not be mixed: 1200
% tonnes of Benzol, 700 tonnes of Butanol, 1000 tonnes of Propanol,
% 450 tonnes of Styrene, and 1200 tonnes of THF. Nine tanks of
% different capacities are available on site. Some of them are
% already partially filled with a liquid. The following table lists
% the characteristics of the tanks (in tonnes). Into which tanks
% should the ships be unloaded (question 1) to maximize the capacity
% of the tanks that remain unused, or (question 2) to maximize the
% number of tanks that remain free?
%
% Characteristics of tanks
% 
% +---------------+---+------+---+---+---+---+---+---+---+
% |Tank           |  1|    2 |  3|  4|  5|  6|  7|  8|  9|
% +---------------+---+------+---+---+---+---+---+---+---+
% |Capacity       |500|  400 |400|600|600|900|800|800|800|
% |Current product|  -|Benzol|  -|  -|  -|  -|THF|  -|  -|
% |Quantity       |  0|  100 |  0|  0|  0|  0|300|  0|  0|
% +---------------+---+------+---+---+---+---+---+---+---+
%
% VARIABLES
%
% capacity                   Capacity of Empty tanks
% products                   Amount of incoming chemicals
%
% RESULTS
%
% In order to interpret the results, try the following:
% Result = tankloadingEx(2);

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
% Written Oct 10, 2005.   Last modified Dec 8, 2005.

function Result = tankloadingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

capacity    = [500;0;400;600;600;900;0;800;800];
products    = [1200-300;700;1000;450;1200-500];

Prob = tankloading(capacity, products);
Result = tomRun('cplex', Prob, PriLev);
Result.remaining = sum(capacity)-Result.f_k;

if PriLev > 1,
   disp('NB: Tank 2 and tank 7 are filled before the optimization starts.')
   temp   = reshape(Result.x_k,9,5);
   ships  = length(products);
   for ship = 1:ships,
      tank = find(temp(:,ship));
      disp(['Ship number ' num2str(ship) ' unloads in tank(s) ' num2str(tank')])
   end  
end

% MODIFICATION LOG
%
% 051010 med   Created
% 051208 med   Added remaining the return
% 060112 per   Added documentation.
% 060124 per   Interpretation of results upgraded.
% 060125 per   Moved disp to end

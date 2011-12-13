%% Tank Loading
%
%% Problem description
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
%  +---------------+---+------+---+---+---+---+---+---+---+
%  |Tank           |  1|    2 |  3|  4|  5|  6|  7|  8|  9|
%  +---------------+---+------+---+---+---+---+---+---+---+
%  |Capacity       |500|  400 |400|600|600|900|800|800|800|
%  |Current product|  -|Benzol|  -|  -|  -|  -|THF|  -|  -|
%  |Quantity       |  0|  100 |  0|  0|  0|  0|300|  0|  0|
%  +---------------+---+------+---+---+---+---+---+---+---+
%
%% Variables
%
%  capacity                   Capacity of Empty tanks
%  products                   Amount of incoming chemicals
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
capacity = [500;0;400;600;600;900;0;800;800];
products = [1200-300;700;1000;450;1200-500];

l = length(products);
m = length(capacity);

load = tom('load',l,m,'int');

% Binary variables
bnds = {0 <= load <= 1};

% Load constraint, the liquids have to go
con1 = {load*capacity >= products};

% Only one product per tank
con2 = {sum(load,1) <= 1};

% Objective
objective = sum(load*capacity);

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Tank Loading';
sol = ezsolve(objective,constraints,[],options);

remaining = sum(capacity)-subs(objective,sol);

PriLev = 1;
if PriLev > 0
    disp('NB: Tank 2 and tank 7 are filled before the optimization starts.')
    temp   = sol.load';
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
% 090308 med   Converted to tomSym
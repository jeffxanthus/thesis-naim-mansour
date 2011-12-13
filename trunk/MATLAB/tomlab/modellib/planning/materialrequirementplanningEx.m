% function Result = materialrequirementplanningEx(PriLev)
%
% Creates a TOMLAB MIP problem for material requirement planning
%
% MATERIAL REQUIREMENT AND PLANNING
%
% The company Minorette produces two types of large toy lorries for
% children: blue removal vans and red tank lorries. Each of these 
% lorries is assembled from thirteen items. See below for the
% breakdown of components and the table below for the prices of the
% components.
%
% Prices of components
%
%+----------+--------------+--------+----------+---------+-----------+
%|  Wheel   |Steel bar     |Bumper  |  Chassis |  Cabin  |Door window|
%+----------+--------------+--------+----------+---------+-----------+
%|  $ 0.30  |   $ 1        |$ 0.20  |  $ 0.80  |  $ 2.75 |  $ 0.10   |
%+----------+--------------+--------+----------+---------+-----------+
%|Windscreen|Blue container|Red tank|Blue motor|Red motor| Headlight |
%+----------+--------------+--------+----------+---------+-----------+
%|  $ 0.29  | $ 2.60       | $ 3    |  $ 1.65  |  $ 1.65 |  $ 0.15   |
%+----------+--------------+--------+----------+---------+-----------+
%
%
% Breakdown of components (Gozinto graph)
%
%                           Blue lorry
%                           (or red)
%                               |
%                               |
%                               |
%   +-----------+---------------+----------+-----------+
%   |           |               |          |           |
%  1|          1|              1|         1|          2|
%   |           |               |          |           |
%  Assembled   Blue container  Assembled  Blue motor  Headlight
%  chassis     (or red tank)   cabin      (or red)  
%   |                             |
%   |                             |
%   +-------+-----+        +------+-------+
%   |       |     |        |      |       |
%  2|      2|    1|       1|     2|      1|
%   |       |     |        |      |       |
%  Bumper  Axle  Chassis  Cabin  Door    Windscreen
%           |                     window
%           |
%        +--+---+
%        |      |
%       2|     1|
%        |      |
%       Wheel  Steel bar
% 
% The subsets (axles, chassis, blue or red cabin) may be assembled by
% the company using the components, or subcontracted. The following 
% table lists the costs for assembling and for subcontracting these
% subsets, together with the capacities of the company. The assembly
% costs do not take into account the buying prices of their
% components.
%
% Subcontracting and assembly costs, assembly capacities
%
% +--------------+------+-----------------+---------------+----------+---------+
% |              | Axle |Assembled chassis|Assembled cabin|Blue lorry|Red lorry|
% +--------------+------+-----------------+---------------+----------+---------+
% |Subcontracting|$12.75|      $ 30       |    $ 3        |     -    |     -   |
% |Assembly      |$6.80 |      $ 3.55     |    $ 3.20     |   $ 2.20 |   $ 2.60|
% |Capacity      | 600  |       4000      |     3000      |    4000  |    5000 |
% +--------------+------+-----------------+---------------+----------+---------+
%
% For the next month, Minorette has the same demand forecast of 3000
% pieces for the two types of lorries. At present, the company has no
% stock. Which quantities of the different items should Minorette buy 
% or subcontract to satisfy the demand while minimizing the 
% production cost?
%
% VARIABLES
%
% demand                     The demand for tank and container lorries 
% compprices                 The price of the components 
% subcontr                   Cost for using a subcontracter to assemble 
% assembly                   Price for own assembly 
% finalassembly              Price for final assembly 
% capacity                   Capacity to assemble components 
% finalcapacity              Capacity to assemble final step 
% assemmat                   12 first columns: component prices
%                            next 3 columns: subcontracter price
%                            next 3 columns: own assembly price
%                            last 2 columns: final assembly
%
% RESULTS
%
% For an interpretation of the results, try:
% Result            = materialrequirementplanningEx(2);
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
% Written Oct 18, 2005.   Last modified Oct 18, 2005.

function Result = materialrequirementplanningEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand          = [3000;3000];
compprices      = [.3;1;.2;.8;2.75;.1;.29;2.6;3;1.65;1.65;.15];
subcontr        = [12.75;30;3];
assembly        = [6.8;3.55;3.20];
finalassembly   = [2.2;2.6];
capacity        = [600;4000;3000];
finalcapacity   = [4000;5000];

%Buy preprod (12), Buy subcontr (3), assemble (3), finalass (2)
assemmat        = [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2  0  0  0  0;...
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  0;...
      0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -2  0  0  0;...
      0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -2  0  0  0;...
      0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0;...
      0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0;...
      0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0 -2  0  0;...
      0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0 -1  0  0;...
      0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1  0 -1 -1;...
      0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0 -1  0;...
      0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0 -1;...
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -1 -1;...
      0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0 -1  0;...
      0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0 -1;...
      0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 -2 -2];

Prob = materialrequirementplanning(demand, compprices, assemmat...
   , subcontr, assembly, finalassembly, capacity, finalcapacity);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   l                 = length(compprices);
   l2                = length(subcontr);
   l3                = length(assembly);
   l4                = length(finalassembly);
   buycomponent      = Result.x_k(1         : l);
   buysubcontracter  = Result.x_k(1+l       : l+l2);
   buyassembly       = Result.x_k(1+l+l2    : l+l2+l3);
   buyfinalassembly  = Result.x_k(1+l+l2+l3 : l+l2+l3+l4);
   disp('Buy these components:')   
   for i = 1:l,
      if buycomponent(i) > 0,
         disp(['   ' num2str(buycomponent(i)) ' units of type ' num2str(i) ])
      end
   end
   disp('Use these subcontracters:')   
   for i = 1:l2,
      if buysubcontracter(i) > 0,
         disp(['   ' num2str(buysubcontracter(i)) ' assemblies of type ' num2str(i) ])
      end
   end
   disp('Do this assembly:')   
   for i = 1:l3,
      if buyassembly(i) > 0,
         disp(['   ' num2str(buyassembly(i)) ' assemblies of type ' num2str(i) ])
      end
   end
   disp('Do this final assembly:')   
   for i = 1:l4,
      if buyfinalassembly(i) > 0,
         disp(['   ' num2str(buyfinalassembly(i)) ' assemblies of type ' num2str(i) ])
      end
   end
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end

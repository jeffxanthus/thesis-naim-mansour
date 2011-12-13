% function Result = planningtheprodofeleccompEx(PriLev)
%
% Creates a TOMLAB MIP problem for planning the production of 
% electronic components
%
% PLANNING THE PRODUCTION OF ELECTRONIC COMPONENTS
%
% To augment its competitiveness a small business wishes to improve
% the production of its best selling products. One of its main
% activities is the production of cards with microchips and 
% electronic badges. The company also produces the components for
% these cards and badges. Good planning of the production of these
% components therefore constitutes a decisive success factor for the
% company. The demand for components is internal in this case and 
% hence easy to anticipate. 
%
% For the next six months the production of four products with
% references X43-M1, X43-M2, Y54-N1, Y54-N2 is to be planned. The 
% production of these components is sensitive to variations of the
% level of production, and every change leads to a non-negligible 
% cost through controls and readjustments. The company therefore
% wishes to minimize the cost associated with these changes whilst
% also taking into account the production and storage costs.
%
% The demand data per time period, the production and storage costs,
% and the initial and desired final stock levels for every product
% are listed in the following table. When the production level
% changes, readjustments of the machines and controls have to be
% carried out for the current month. The cost incurred is 
% proportional to the increase or reduction of the production
% compared to the preceding month. The cost for an increase of the 
% production is $ 1 per unit but only $ 0.50 for a decrease of the
% production level.
%
% Data for the four products
%
% +------------------------------------+------------------+-------------+
% |           Demands                  |        Cost      |     Stock   |
% +------+----+----+----+----+----+----+----------+-------+-------+-----+
% | Month|  1 |  2 |  3 |  4 |  5 |  6 |Production|Storage|Initial|Final|
% +------+----+----+----+----+----+----+----------+-------+-------+-----+
% |X43-M1|1500|3000|2000|4000|2000|2500|    20    |  0.4  |  10   | 50  |
% |X43-M2|1300| 800| 800|1000|1100| 900|    25    |  0.5  |   0   | 10  |
% |Y54-N1|2200|1500|2900|1800|1200|2100|    10    |  0.3  |   0   | 10  |
% |Y54-N2|1400|1600|1500|1000|1100|1200|    15    |  0.3  |   0   | 10  |
% +------+----+----+----+----+----+----+----------+-------+-------+-----+
%
% What is the production plan that minimizes the sum of costs
% incurred through changes of the production level, production
% and storage costs?
%
% VARIABLES
%
% demand                      the demand for each component and month
% prodcost                    Cost to produce a component 
% storagecost                 Cost to store a component 
% initialstock                Initial stock 
% finalstock                  Final stock 
% increasecost, decreasecost  Increase or decrease in cost
%                             when producing more or less this month
%                             than last month.
%
% RESULTS
%
% for an interpretation of the results, try:
% Result      = planningtheprodofeleccompEx(2);
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
% Written Oct 18, 2005.   Last modified Feb 3, 2006.

function Result = planningtheprodofeleccompEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand          = [1500 3000 2000 4000 2000 2500;...
      1300  800  800 1000 1100  900;...
      2200 1500 2900 1800 1200 2100;...
      1400 1600 1500 1000 1100 1200];
prodcost        = [20;25;10;15];
storagecost     = [.4;.5;.3;.3];
initialstock    = [10;0;50;0];
finalstock      = [50;10;30;10];
increasecost    = 1;
decreasecost    = 0.5;

Prob = planningtheprodofeleccomp(demand, prodcost,...
   storagecost, initialstock, finalstock, increasecost, decreasecost);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   l1          = 10; % 4 components + 4 stockslots + increase + decrease 
   l2          =  6; % the number of time segments
   temp        = reshape(Result.x_k,l1,l2);
   for month = 1:l2,
      disp(['Solution for month ' num2str(month)]) 
      disp(['produce  ' num2str(temp( 1: 4,month)')])
      disp(['stock    ' num2str(temp( 5: 8,month)')])
      disp(['increase ' num2str(temp( 9: 9,month)')])
      disp(['decrease ' num2str(temp(10:10,month)')])
      disp(' ')
   end
   
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
% 060203 med   Removed printing of temp
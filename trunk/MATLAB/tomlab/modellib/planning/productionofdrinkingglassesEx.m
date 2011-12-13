% function Result = productionofdrinkingglassesEx(PriLev)
%
% Creates a TOMLAB MIP problem for production of drinking glasses
%
% PRODUCTION OF DRINKING GLASSES
%
% The main activity of a company in northern France is the production
% of drinking glasses. It currently sells six different types
% (V1 to V6), that are produced in batches of 1000 glasses, and wishes
% to plan its production for the next 12 weeks. The batches may be 
% incomplete (fewer than 1000 glasses). The demand in thousands for
% the 12 coming weeks and for every glass type is given in the 
% following table.
%
% Demands for the planning period (batches of 1000 glasses)
%
% +----+--+--+--+--+--+--+--+--+--+--+--+--+
% |Week| 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|
% +----+--+--+--+--+--+--+--+--+--+--+--+--+
% |V1  |20|22|18|35|17|19|23|20|29|30|28|32|
% |V2  |17|19|23|20|11|10|12|34|21|23|30|12|
% |V3  |18|35|17|10| 9|21|23|15|10| 0|13|17|
% |V4  |31|45|24|38|41|20|19|37|28|12|30|37|
% |V5  |23|20|23|15|10|22|18|30|28| 7|15|10|
% |V6  |22|18|20|19|18|35| 0|28|12|30|21|23|
% +----+--+--+--+--+--+--+--+--+--+--+--+--+
%
% For every glass type the initial stock is known, as well as the
% required final stock level (in thousands). Per batch of every glass
% type, the production and storage costs in $ are given, together
% with the required working time for workers and machines (in hours),
% and the required storage space (measured in numbers of trays).
%
% The number of working hours of the personnel is limited to 450  
% hours per week, and the machines have a weekly capacity of 850
% hours. Storage space for up to 1000 trays is available. Which
% quantities of the different glass types need to be produced in 
% every period to minimize the total cost of production and storage?
%
% Data for the six glass types
%
% +--+----------+-------+-------+-----+----------+-----------+-------+
% |  |Production|Storage|Initial|Final|          |           |Storage|
% |  |  cost    |  cost | stock |stock|Timeworker|Timemachine| space |
% +--+----------+-------+-------+-----+----------+-----------+-------+
% |V1|   100    |   25  |   50  | 10  |    3     |    2      |   4   |
% |V2|    80    |   28  |   20  | 10  |    3     |    1      |   5   |
% |V3|   110    |   25  |    0  | 10  |    3     |    4      |   5   |
% |V4|    90    |   27  |   15  | 10  |    2     |    8      |   6   |
% |V5|   200    |   10  |    0  | 10  |    4     |   11      |   4   |
% |V6|   140    |   20  |   10  | 10  |    4     |    9      |   9   |
% +--+----------+-------+-------+-----+----------+-----------+-------+
%
% VARIABLES
%
% demand                     Weekly demand of V1-V6
% workermax                  The staffs weekly capacity in hours
% machinemax                 The machines weekly capacity in hours
% maxstorage                 Maximal storage space
% productioncost             Cost to produce a batch of a glasstype.
% storagecost                Cost to store a batch of a glasstype
% initialstock               Initial stock of glasstypes
% finalstock                 Required final stock
% timeworker                 Personnel-time required for a batch
% timemachine                Machine-time required for a batch
% storagespace               Space required by batch in storage
%
% RESULTS
%
% for an interpretation of the results, try
% Result = productionofdrinkingglassesEx(2);
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
% Written Oct 17, 2005.   Last modified Oct 17, 2005.

function Result = productionofdrinkingglassesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand          = [20 22 18 35 17 19 23 20 29 30 28 32;...
      17 19 23 20 11 10 12 34 21 23 30 12;...
      18 35 17 10  9 21 23 15 10  0 13 17;...
      31 45 24 38 41 20 19 37 28 12 30 37;...
      23 20 23 15 10 22 18 30 28  7 15 10;...
      22 18 20 19 18 35  0 28 12 30 21 23];

workermax       = 450; % Modified from 390, otherwise infeasible
machinemax      = 850;
maxstorage      = 1000;

productioncost  = [100;80;110;90;200;140];
storagecost     = [25;28;25;27;10;20];
initialstock    = [50;20;0;15;0;10];
finalstock      = [10;10;10;10;10;10];
timeworker      = [3;3;3;2;4;4];
timemachine     = [2;1;4;8;11;9];
storagespace    = [4;5;5;6;4;9];

Prob = productionofdrinkingglasses(demand, workermax, machinemax, maxstorage...
   , productioncost, storagecost, initialstock, finalstock, timeworker...
   , timemachine, storagespace);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   
   [g,w] = size(demand); % g = glass_types, W = weeks
   
   temp = reshape(Result.x_k,g,2,w);
   
   for i = 1:w,
      disp(['results for week ' num2str(i) ':'])
      for j = 1:g,
         if temp(j,1,i) > 0,
            disp(['   ' num2str(temp(j,1,i)) ' batches of type ' num2str(j) ' should be produced' ])
         end
         if temp(j,2,i) > 0,
            disp(['   ' num2str(temp(j,2,i)) ' batches of type ' num2str(j) ' should be stored to next month' ])
         end
      end
   end
end

% MODIFICATION LOG
%
% 051017 med   Created.
% 060109 per   Added documentation.
% 060125 per   Moved disp to end

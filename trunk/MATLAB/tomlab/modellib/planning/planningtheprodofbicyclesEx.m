% function Result = planningtheprodofbicyclesEx(PriLev)
%
% Creates a TOMLAB MIP problem for planning the production of bicycles
%
% PLANNING THE PRODUCTION OF BICYCLES
% 
% A company produces bicycles for children. The sales forecast in
% thousand of units for the coming year are given in the following
% table. The company has a capacity of 30,000 bicycles per month.
% It is possible to augment the production by up to 50% through
% overtime working, but this increases the production cost for a
% bicycle from the usual $ 32 to $ 40. 
%
% Sales forecasts for the coming year in thousand units
% 
% +---+---+---+---+---+---+---+---+---+---+---+---+
% |Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|
% +---+---+---+---+---+---+---+---+---+---+---+---+
% | 30| 15| 15| 25| 33| 40| 45| 45| 26| 14| 25| 30|
% +---+---+---+---+---+---+---+---+---+---+---+---+
%
% Currently there are 2,000 bicycles in stock. The storage costs have
% been calculated as $ 5 per unit held in stock at the end of a month.
% We assume that the storage capacity at the company is virtually
% unlimited (in practice this means that the real capacity, that is 
% quite obviously limited, does not impose any limits in our case).
% We are at the first of January. Which quantities need to be
% produced and stored in the course of the next twelve months in order
% to satisfy the forecast demand and minimize the total cost?
%
% VARIABLES
%
% normcapacity               the normal production capacity
% extracapacity              extra capacity
% normcost                   normal cost
% extracost                  cost per bike if overtime
% demand                     bikes wanted per month
% startstock                 bikes in store 
% storagecost                cost to have a bike in store one month
%
% RESULTS
%
% For an interpretation of the results, try this:
% Result = planningtheprodofbicyclesEx(2);
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

function Result = planningtheprodofbicyclesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

normcapacity    = 30000;
extracapacity   = 15000;
normcost        = 32;
extracost       = 40;
demand          = [30;15;15;25;33;40;45;45;26;14;25;30]*1000;
startstock      = 2000;
storagecost     = 5;

Prob = planningtheprodofbicycles(normcapacity, extracapacity, normcost,...
   extracost, demand, startstock, storagecost);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   months = length(demand);
   temp = reshape(Result.x_k,months,3);
   
   
   
   for i = 1:months,
      disp(['Month ' num2str(i) ':'])
      disp(['   produce ' num2str(temp(i,1)) ' regular bikes' ])
      if temp(i,2) > 0,
         disp(['   and     ' num2str(temp(i,2)) ' extras' ])
      end
      if temp(i,3) > 0,
         disp(['   let     ' num2str(temp(i,3)) ' be stored' ])
      end
   end
  
end

% MODIFICATION LOG
%
% 051017 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end

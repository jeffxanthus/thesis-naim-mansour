% function Result = productionofelectricityEx(PriLev)
%
% Creates a TOMLAB MIP problem for production of electricity
%
% PRODUCTION OF ELECTRICITY
%
% Power generators of four different types are available to satisfy the
% daily electricity demands (in megawatts) summarized in the following 
% table. We consider a sliding time horizon: the period 10pm-12am of day d
% is followed by the period 0am-6am of day d + 1. 
%
% Daily electricity demands (in MW)
%
% +-------+-------+-------+--------+--------+-------+--------+---------+
% |Period |0am-6am|6am-9am|9am-12pm|12pm-2pm|2pm-6pm|6pm-10pm|10pm-12am|
% +-------+-------+-------+--------+--------+-------+--------+---------+
% |Demand |  12000|  32000|   25000|   36000|  25000|   30000|    18000|
% +-------+-------+-------+--------+--------+-------+--------+---------+
%
% The power generators of the same type have a maximum capacity and may be
% connected to the network starting from a certain minimal power output. 
% They have a start-up cost, a fixed hourly cost for working at minimal
% power, and an hourly cost per additional megawatt for anything beyond
% the minimal output. These data are given in the following table. 
%
% Description of power generators
% +---------+-----------+------+--------+--------+------------+--------+
% |Available|Min. output|Max.  |capacity|Fix cost|Add. MW cost|Start-up|
% |number   |in MW      |in MW |$/h     |        |$/h         |cost    |
% +---------+-----------+------+--------+--------+------------+--------+
% |Type 1   |  10       |  750 |  1750  |  2250  |  2.7       | 5000   |
% |Type 2   |   4       | 1000 |  1500  |  1800  |  2.2       | 1600   |
% |Type 3   |   8       | 1200 |  2000  |  3750  |  1.8       | 2400   |
% |Type 4   |   3       | 1800 |  3500  |  4800  |  3.8       | 1200   |
% +---------+-----------+------+--------+--------+------------+--------+
% 
% A power generator can only be started or stopped at the beginning of a
% time period. As opposed to the start, stopping a power plant does not
% cost anything. At any moment, the working power generators must be able
% to cope with an increase by 20% of the demand forecast. Which power
% generators should be used in every period in order to minimize the total
% daily cost?
%
% VARIABLES
%
% demand                     MW needed each period
% available                  Generators of each type available
% mincap                     Minimal output of generator
% maxcap                     Maximal output of generator
% fixcost                    Fix cost per hour
% runningcost                Cost per hour and MW above mincap
% startcost                  Cost for starting a generator
% periodlengths              Hours in each period
%
% RESULTS
%
% For an interpretation of the results use PriLev > 1, for example:
% Result = productionofelectricityEx(2);
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
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Result = productionofelectricityEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand      = [12;32;25;36;25;30;18]*1000;
available   = [10;4;8;3];
mincap      = [750;1000;1200;1800];
maxcap      = [1750;1500;2000;3500];
fixcost     = [2250;1800;3750;4800];
runningcost = [2.7;2.2;1.8;3.8];
startcost   = [5000;1600;2400;1200];
periodlengths = [6;3;3;2;4;4;2];
Prob = productionofelectricity(demand, available, mincap, maxcap,...
   fixcost, runningcost, startcost, periodlengths);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   temp      = reshape(Result.x_k,length(available(:)),length(demand(:)),3);
   intervals = ['00-06';'06-09';'09-12';'12-14';'14-18';'18-22';'22-24'];
   starts    = temp(:,:,1);
   running   = temp(:,:,2);
   extra     = temp(:,:,3);
   [reacs, times] = size(starts);
   for time = 1:times,
      disp([' At ' intervals(time,:) '...'])
      for reac = 1:reacs,
         if running(reac,time) > 0,
            disp([   '   there are ' num2str(running(reac,time))  ...
                     ' reactors of type ' num2str(reac) ' running' ])
            if starts(reac,time) > 0,
               disp(['      (we started ' num2str( starts(reac,time))...
                     ' additional reactors)'])
            end
            if extra(reac,time) > 0,
               disp(['      (we also produce ' num2str(  extra(reac,time))...
                     ' MW more than the minimal level)'])
            end
         end
      end
   end
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060109 per   Added documentation.
% 060125 per   Moved disp to end

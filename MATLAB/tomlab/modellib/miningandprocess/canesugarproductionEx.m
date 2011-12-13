% function Result = canesugarproductionEx(PriLev)
%
% Creates a TOMLAB LP problem for cane sugar production
%
% CANE SUGAR PRODUCTION
%
% PROBLEM
%
% The harvest of cane sugar in Australia is highly mechanized. The sugar
% cane is immediately transported to a sugar house in wagons that run on
% a network of small rail tracks. The sugar content of a wagon load
% depends on the field it has been harvested from and on the maturity of
% the sugar cane. Once harvested, the sugar content decreases rapidly
% through fermentation and the wagon load will entirely lose its value
% after a certain time. At this moment, eleven wagons all loaded with the
% same quantity have arrived at the sugar house. They have been examined
% to find out the hourly loss and the remaining life span (in hours) of 
% every wagon, these data are summarized in the following table.
%
% Table: Properties of the lots of cane sugar
%
% Lot             1  2  3  4  5  6  7  8  9 10 11
% Loss (kg/h)    43 26 37 28 13 54 62 49 19 28 30
% Life span (h)   8  8  2  8  4  8  8  8  8  8  8
%
% Every lot may be processed by any of the three, fully equivalent
% production lines of the sugar house. The processing of a lot takes two
% hours. It must be finished at the latest at the end of the life span of
% the wagon load. The manager of the sugar house wishes to determine a 
% production schedule for the currently available lots that minimizes the
% total loss of sugar.
%
% VARIABLES
%
% proclines                  the number of processing lines
% proctime                   time in hours to process a wagon
% loss                       loss of sugar in kg per hour
% lifespan                   how long the sugar in a wagon will last
%
% RESULTS
%
% For an interpretation of the results, set PriLev > 1, for example:
% Result = canesugarproductionEx(2);
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
%   Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Dec 8, 2005.

function Result = canesugarproductionEx(PriLev)

if nargin < 1
   PriLev = 1;
end

proclines   = [3];
proctime    = [2];
loss        = [43;26;37;28;13;54;62;49;19;28;30];
lifespan    = [8;8;2;8;4;8;8;8;8;8;8];

Prob = canesugarproduction(proclines, proctime, loss, lifespan);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,   
   lots      = length(lifespan);
   timeslots = ceil(lots/proclines(1));  % even a fraction means work
   temp      = reshape(Result.x_k,lots,timeslots);
   
   disp(['To minimize loss (' num2str(Result.f_k) ') use this schema'])
   
   for time = 1:timeslots,
      disp(['at ' num2str((time-1)*proctime(1)+8) '.00:'])
      idx = find(temp(:,time));
      disp(['   process the lots ' num2str(idx') ])
   end
   disp(['at ' num2str((time)*proctime(1)+8) '.00:'])
   disp( '   harvesting completed')   
end

% MODIFICATION LOG
%
% 051007 med   Created
% 051208 med   Lifespan factor wrong
% 060109 per   lifespan(9) changed to 8
% 060109 per   Added documentation.
% 060125 per   Moved disp to end

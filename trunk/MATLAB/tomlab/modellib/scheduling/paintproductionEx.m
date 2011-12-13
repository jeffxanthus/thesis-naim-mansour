% function Result = paintproductionEx(PriLev)
%
% Creates a TOMLAB MIP problem for paint production
%
%
% PAINT PRODUCTION
%
% As a part of its weekly production a paint company produces five
% batches of paints, always the same, for some big clients who have a
% stable demand. Every paint batch is produced in a single production
% process, all in the same blender that needs to be cleaned between
% two batches. The durations of blending paint batches 1 to 5 are
% respectively 40, 35, 45, 32, and 50 minutes. The cleaning times
% depend on the colors and the paint types. For example, a long
% cleaning period is required if an oil-based paint is produced after
% a water-based paint, or to produce white paint after a dark color.
% The times are given in minutes in the following table CLEAN where
% CLEANij denotes the cleaning time between batch i and batch j. 
%
% Matrix of cleaning times
%
% +-+--+--+--+--+--+
% | | 1| 2| 3| 4| 5|
% +-+--+--+--+--+--+
% |1| 0|11| 7|13|11|
% |2| 5| 0|13|15|15|
% |3|13|15| 0|23|11|
% |4| 9|13| 5| 0| 3|
% |5| 3| 7| 7| 7| 0|
% +-+--+--+--+--+--+
%
% Since the company also has other activities, it wishes to deal
% with this weekly production in the shortest possible time (blending
% and cleaning). Which is the corresponding order of paint batches?
% The order will be applied every week, so the cleaning time between
% the last batch of one week and the first of the following week
% needs to be counted for the total duration of cleaning.
%
% VARIABLES
%
% cleantimes                 Times to clean from batch i to j
% prodtimes                  Production times per batch
%
% RESULTS
%
% For an interpretation of the results, use PriLev > 1, for example:
% Result  = paintproductionEx(2);       
%
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
% Written Oct 17, 2005.   Last modified Jan 2, 2006.

function Result = paintproductionEx(PriLev)

if nargin < 1
   PriLev = 1;
end

cleantimes = [ 0 11  7 13 11;...
      5  0 13 15 15;...
      13 15  0 23 11;...
      9 13  5  0  3;...
      3  7  7  7  0];

prodtimes  = [40;35;45;32;50];            

Prob = paintproduction(prodtimes, cleantimes);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   batches = size(cleantimes,1);                    % number of batches
   temp1   = reshape(Result.x_k,batches,batches+1); % reshape x_K
   link    = [];                                    % connections
   
   for i = 1:batches,
      for j = 1:batches,
         if temp1(i,j) == 1
            link = [[i j ]; link ];                 % finding connections
            
         end
      end
   end
   
   first = link(1:1);                               % start batch
   next =  link(1,2);                               % next batch
   order = [first];                                 % ordered batches
   
   for k = 1:batches,
      order = [order next];                         % adding next
      next =  link(find(link(:,1)==next),2);        % finding new next
   end
   disp(['one best order: ' num2str(order)])        % display solution
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060111 per   Added documentation.
% 060126 per   Moved disp to end

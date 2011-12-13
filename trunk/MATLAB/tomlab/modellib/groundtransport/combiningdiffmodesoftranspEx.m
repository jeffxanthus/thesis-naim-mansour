% function Result = combiningdiffmodesoftranspEx(PriLev)
%
% Creates a TOMLAB MIP problem for combining different modes of transport
%
% COMBINING DIFFERENT MODES OF TRANSPORT
%
% A load of 20 tonnes needs to be transported on a route passing
% through five cities, with a choice of three different modes of
% transport: rail, road, and air. In any of the three intermediate
% cities it is possible to change the mode of transport but the load
% uses a single mode of transport between two consecutive cities.
% The following table lists the cost of transport in $ per tonne
% between the pairs of cities.
%
% Transport costs with different modes
%
%
%  Pairs of cities
% +----+---+---+---+---+
% |    |1–2|2–3|3–4|4–5|
% +----+---+---+---+---+
% |Rail| 30| 25| 40| 60|
% |Road| 25| 40| 45| 50|
% |Air | 40| 20| 50| 45|
% +----+---+---+---+---+
%
% The next table summarizes the costs for changing the mode of
% transport in $ per tonne. The cost is independent of location.
%
% 
% Cost for changing the mode of transport
%
% +---------+----+----+---+
% |from \ to|Rail|Road|Air|
% +---------+----+----+---+
% |Rail     |  0 |  5 | 12|
% |Road     |  8 |  0 | 10|
% |Air      | 15 | 10 |  0|
% +---------+----+----+---+
%
% How should we organize the transport of the load at the least cost?
%
% VARIABLES
% 
% transpcost                 Transport costs
% changecost                 Cost to change mode of transport
% demand                     Load to transport 
%
% RESULTS
%
% For an interpretation of the results, try this:
% Result = combiningdiffmodesoftranspEx(2);
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
% Written Oct 20, 2005.   Last modified Oct 20, 2005.

function Result = combiningdiffmodesoftranspEx(PriLev)

if nargin < 1
   PriLev = 1;
end

transpcost     = [ 30 25 40 60 ;...
      25 40 45 50 ;...
      40 20 50 45];

changecost     = [  0  5  12;...
      8  0  10;...
      15 10   0];

demand         = 30;                

Prob = combiningdiffmodesoftransp(transpcost, changecost, demand);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   modes  = reshape(Result.x_k(1:3*4),3,4);
   means  = ['rail';'road';'air '];
   disp(['the min cost (' num2str(Result.f_k) ') is found by,'])
   for m = 1:size(modes,2),
      disp(['   going from town ' num2str(m) ' to town ' num2str(m+1) ...
            ' by ' means(find(modes(:,m)),:)])
   end
end

% MODIFICATION LOG
%
% 051020 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end

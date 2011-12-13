% function Result = productionofalloysEx(PriLev)
%
% Creates a TOMLAB LP problem for production design of alloys
%
% PRODUCTION OF ALLOYS
%
% The company Steel has received an order for 500 tonnes of steel to
% be used in shipbuilding. This steel must have the following 
% characteristics (‘grades’). 
%
% Characteristics of steel ordered
% +---------------+-------+-------+
% |Chemical       |Minimum|Maximum|
% |Element        |Grade  |Grade  |
% +---------------+-------+-------+
% |Carbon (C)     |  2.0  | 3.0   |
% |Copper (Cu)    |  0.4  | 0.6   |
% |Manganese (Mn) |  1.2  | 1.65  |
% +---------------+-------+-------+
%
% The company has seven different raw materials in stock that may be
% used for the production of this steel. The table below lists the
% grades, available amounts and prices for all raw materials.
%
% Raw material grades, availabilities, and prices
% +----------------+---+----+----+------------+----+
% |Raw material    |C %|Cu %|Mn %|Availability|Cost|
% +----------------+---+----+----+------------+----+
% |Iron alloy 1    |2.5| 0  |1.3 |     400    |200 |
% |Iron alloy 2    | 3 | 0  |0.8 |     300    |250 |
% |Iron alloy 3    | 0 | 0.3|0   |     600    |150 |
% |Copper alloy 1  | 0 | 90 |0   |     500    |220 |
% |Copper alloy 2  | 0 | 96 |4   |     200    |240 |
% |Aluminum alloy 1| 0 | 0.4|1.2 |     300    |200 |
% |Aluminum alloy 2| 0 | 0.6|0   |     250    |165 |
% +----------------+---+----+----+------------+----+
%
% The objective is to determine the composition of the steel that
% minimizes the production cost.
%
% VARIABLES
% 
% compsize                   Amount of steel to produce
% mincomp                    Minimal contents of each component.
% maxcomp                    Maximal contents of each component.
% rawcompmat                 Contents of raw material.
% rawavail                   Amount available of raw material.
% rawcost                    Cost of raw material.
%
% RESULTS
%
% Result.x_k has the same length as there are raw material. In it we
% find how much of each component we should buy. For an 
% interpretation of the results, use PriLev > 1, for example:
% Result = productionofalloysEx(2);
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
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Result = productionofalloysEx(PriLev)

if nargin < 1
   PriLev = 1;
end

compsize   = 500;
mincomp    = [2;.4;1.2];
maxcomp    = [3;.6;1.65];
rawcompmat = [[2.5;  3;  0;  0;  0;  0;  0],...
      [  0;  0; .3; 90; 96; .4; .6],...
      [1.3; .8;  0;  0;  4;1.2;  0]];
rawavail   = [ 400;300;600;500;200;300;250];
rawcost    = [ 200;250;150;220;240;200;165];

Prob = productionofalloys(compsize, mincomp, maxcomp, rawcompmat, rawavail, rawcost);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   names = [ 'Iron 1' ; 'Iron 2' ; 'Iron 3' ; 'Copp 1' ; ...
         'Copp 2' ; 'Alum 1' ; 'Alum 2']; 
   disp(['the optimal steel contains:'])
   for i = 1:length(Result.x_k),
      if Result.x_k(i) > 0,
         disp(['   ' num2str(Result.x_k(i)) ' tonnes of ' names(i,:) ])
      end
   end
end

% MODIFICATION LOG
%
% 051007 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end

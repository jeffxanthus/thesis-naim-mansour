% function Result = carrentalEx(PriLev)
%
% Creates a TOMLAB MIP problem for car rental
%
% CAR RENTAL
%
% A small car rental company has a fleet of 94 vehicles distributed
% among its 10 agencies. The location of every agency is given by its
% geographical coordinates X and Y in a grid based on kilometers. We
% assume that the road distance between agencies is approximately 1.3
% times the Euclidean distance (as the crow flies). The following
% table indicates the coordinates of all agencies, the number of cars
% required the next morning, and the stock of cars in the evening
% preceding this day.
% 
% Description of the vehicle rental agencies
%
% +-------------+--+--+--+--+--+--+--+--+--+--+
% |Agency       | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|
% +-------------+--+--+--+--+--+--+--+--+--+--+
% |X coordinate | 0|20|18|30|35|33| 5| 5|11| 2|
% |Y coordinate | 0|20|10|12| 0|25|27|10| 0|15|
% +-------------+--+--+--+--+--+--+--+--+--+--+
% |Required cars|10| 6| 8|11| 9| 7|15| 7| 9|12|
% |Cars present | 8|13| 4| 8|12| 2|14|11|15| 7|
% +-------------+--+--+--+--+--+--+--+--+--+--+
%               
% Supposing the cost for transporting a car is $ 0.50 per km,
% determine the movements of cars that allow the company to 
% re-establish the required numbers of cars at all agencies,
% minimizing the total cost incurred for transport.
%
% VARIABLES
%
% demand                     Required cars per agency
% stock                      Cars present
% cost                       Cost per km to transport a car
% xcord                      The X-coordinate of agencies
% ycord                      The Y-coordinate of agencies
% n                          Number of agencies
% distance                   A matrix of distances
%
% RESULTS
%
% For an interpretation of the results, try this:
% Result          = carrentalEx(2);
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
% Written Oct 19, 2005.   Last modified Oct 19, 2005.

function Result = carrentalEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand          = [10   6   8  11   9   7  15   7   9  12]';
stock           = [ 8  13   4   8  12   2  14  11  15   7]';
cost            = 0.50;

xcord           = [ 0  20  18  30  35  33   5   5  11   2]';
ycord           = [ 0  20  10  12   0  25  27  10   0  15]';
n               = length(xcord);

distance        = zeros(n,n);

for i=1:n
   for j=1:n
      distance(i,j) = 1.3*sqrt( (xcord(i) - xcord(j))^2 + (ycord(i) - ycord(j))^2);
   end
end

Prob = carrental(cost, distance, stock, demand);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   idx_excess      = find(stock-demand > 0);
   idx_need        = find(stock-demand < 0);
   n_excess        = length(idx_excess);
   n_need          = length(idx_need);
   temp            = reshape(Result.x_k,n_excess,n_need);
   
   disp('THE SENDING OF CARS')
   for i = 1:n_excess,       % scan all positions, disp interpretation
      disp(['agency ' num2str(idx_excess(i)) ' sends: ' ])
      for j = 1:n_need,
         if temp(i,j) ~= 0
            disp(['   ' num2str(temp(i,j)) ' car(s) to agency ' ...
                  num2str(idx_need(j))])
         end
      end
      disp(' ')
   end
   
   disp('THE GETTING OF CARS')
   for j = 1:n_need,
      disp(['agency ' num2str(idx_need(j)) ' gets: ' ])
      for i = 1:n_excess,       % scan all positions, disp interpretation
         if temp(i,j) ~= 0
            disp(['   ' num2str(temp(i,j)) ' car(s) from agency ' ...
                  num2str(idx_excess(i))])
         end
      end
      disp(' ')
   end
end
% MODIFICATION LOG
%
% 051019 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end

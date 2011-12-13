% function Result = planningtheprodoffiberglassEx(PriLev)
%
% Creates a TOMLAB MIP problem for planning the production of 
% fiber glass
%
% PLANNING THE PRODUCTION OF FIBERGLASS
%
% A company produces fiberglass by the cubic meter and wishes to plan
% its production for the next six weeks. The production capacity is
% limited, and this limit takes a different value in every time
% period. The weekly demand is known for the entire planning period.
% The production and storage costs also take different values
% depending on the time period. All data are listed in the following
% table.
%
% Data per week
%
% +----+------------+------+----------+-----------+
% |    | Production |Demand|Production| Storage   |
% |Week|capacity(m3)| (m3) |cost($/m3)|cost ($/m3)|
% +----+------------+------+----------+-----------+
% |  1 |    140     | 100  |   5      |   0.2     |
% |  2 |    100     | 120  |   8      |   0.3     |
% |  3 |    110     | 100  |   6      |   0.2     |
% |  4 |    100     |  90  |   6      |   0.25    |
% |  5 |    120     | 120  |   7      |   0.3     |
% |  6 |    100     | 110  |   6      |   0.4     |
% +----+------------+------+----------+-----------+
% 
% Which is the production plan that minimizes the total cost of
% production and storage?
%
% VARIABLES
%
% capacity                   Production capacity over time
% demand                     Demand over time
% prodcost                   Cost to produce over time
% storcost                   Cost to store over time
%
% RESULTS
%
% For an interpretation of the results, try:
% Result = planningtheprodoffiberglassEx(2);
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

function Result = planningtheprodoffiberglassEx(PriLev)

if nargin < 1
   PriLev = 1;
end

capacity        = [140;100;110;100;120;100];
demand          = [100;120;100; 90;120;110];
prodcost        = [  5;  8;  6;  6;  7;  6];
storcost        = [ .2; .3; .2;.25; .3; .4];

Prob = planningtheprodoffiberglass(capacity, demand,...
   prodcost, storcost);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   times = length(capacity);
   
   produce = Result.x_k(1:times);
   store   = [Result.x_k(times+1:end)' 0]';
   for t = 1:times,
      if produce(t) > 0 | store(t) > 0,
         disp(['during month ' num2str(t) ])
         if produce(t) > 0,
            disp(['   produce  ' num2str(produce(t))])
         end
         if store(t) > 0,
            disp(['   store  ' num2str(store(t)) ' to next month'])
         end
      end
   end
end

% MODIFICATION LOG
%
% 051018 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end

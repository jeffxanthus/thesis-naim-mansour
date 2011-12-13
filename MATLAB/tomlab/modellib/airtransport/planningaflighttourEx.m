% function Result = planningaflighttourEx(PriLev)
%
% Creates a TOMLAB MIP problem for planning a flight tour
%
% PLANNING A FLIGHT TOUR
%
% A country in south-east Asia is experiencing widespread flooding.
% The government, with international help, decides to establish a
% system of supply by air. Unfortunately, only seven runways are
% still in a usable state, among which is the one in the capital. 
%
% The government decides to make the planes leave from the capital,
% have them visit all the other six airports and then come back to
% the capital. The following table lists the distances between the
% airports. Airport A1 is the one in the capital. In which order
% should the airports be visited to minimize the total distance
% covered?
%
% Distance matrix between airports (in km)
%
% +--+---+---+---+---+---+---+
% |  | A2| A3| A4| A5| A6| A7|
% +--+---+---+---+---+---+---+
% |A1|786|549|657|331|559|250|
% |A2|668|979|593|224|905|   |
% |A3|316|607|472|467|   |   |
% |A4|890|769|400|   |   |   |
% |A5|386|559|   |   |   |   |
% |A6|681|   |   |   |   |   |
% +--+---+---+---+---+---+---+
%
% VARIABLES
%
% distance                   The distance matrix
%
% RESULTS
%
% For an interpretation of the results set PriLev > 1, for example:
% Result = planningaflighttourEx(2);
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
% Written Nov 7, 2005.   Last modified Nov 7, 2005.

function Result = planningaflighttourEx(PriLev)

if nargin < 1
   PriLev = 1;
end

distance      = [   0   786   549    657    331    559   250;...
      786     0   668    979    593    224   905;...
      549   668     0    316    607    472   467;...
      657   979   316      0    890    769   400;...
      331   593   607    890      0    386   559;...
      559   224   472    769    386      0   681;...
      250   905   467    400    559    681     0];

Prob = planningaflighttour(distance);

stop = 0;
n = size(distance,1);
while stop<100
   Result = tomRun('cplex', Prob, PriLev);
   % Find a subtour
   idx = zeros(n,n);
   idx(1,1:n) = 1:7;
   for i=1:n
      for j=2:n
         idx(j,i) = find(Result.x_k( (idx(j-1,i)-1)*n+1: idx(j-1,i)*n) ~= 0);
         if idx(j,i) == i
            break
         end
      end
   end
   shortmat = spones(idx);
   len = full(sum(shortmat,1));
   [val, idx2] = min(len);
   idx3 = idx(1:val, idx2);
   
   Aextra = zeros(1,Prob.N);
   for k=1:val-1
      Aextra(1,idx3(k)+n*(idx3(k+1)-1)) = 1;
   end
   
   Prob.A = [Prob.A;Aextra];
   Prob.b_L = [Prob.b_L;-inf];
   Prob.b_U = [Prob.b_U;1];
   Prob.mLin = Prob.mLin+1;
   Result.tour = idx(:,1);
   if val == n
      break
   end
   stop = stop + 1;
end

if PriLev > 1,
   a              =  7;                           % number of airports
   temp           = reshape(Result.x_k,a,a);      % reshape results
   not_home_again = true;                         % help to loop
   current_pos    = 1;                            % we start in 1
   route          = [];                           % initially empty
   route          = [route current_pos];          % add position
   while not_home_again,                          % loop if not home
      current_pos = find(temp(current_pos,:));    % go to next pos
      route = [route current_pos];                % add to route
      if current_pos == 1,                        % if home:
         not_home_again = false;                  %  stop loop
      end
   end
   disp(['Shortest route is: ' num2str(route)])
   disp(['               or: ' num2str(route(end:-1:1))])
   
end

% MODIFICATION LOG
%
% 051107 med   Created.
% 060116 per   Added documentation.
% 060125 per   Moved disp to end

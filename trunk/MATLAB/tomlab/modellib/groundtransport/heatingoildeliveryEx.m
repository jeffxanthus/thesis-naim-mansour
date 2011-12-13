% function Result = heatingoildeliveryEx(PriLev)
%
% Creates a TOMLAB MIP problem for heating oil delivery
%
% HEATING OIL DELIVERY
%
% A transporter has to deliver heating oil from the refinery at
% Donges to a certain number of clients in the west of France. His
% clients are located at Brain-sur-l’Authion, Craquefou, Guérande, la
% Haie Fouassière, Mésanger and les Ponts-de-Cé. The following table
% lists the demands in liters for the different sites. 
%
% Demands by clients (in liters)
% 
% +-----------+-----+-----+-----+-----+
% | BslA| Craq| Guér| H Fo| Mésa| P dC|
% +-----+-----+-----+-----+-----+-----+ 
% |14000| 3000| 6000|16000|15000| 5000|
% +-----+-----+-----+-----+-----+-----+
% 
% The next table contains the distance matrix between the clients and
% the refinery.
%
% Distance matrix (in km)
% +-------------------+----+----+----+----+----+----+----+                                                              
% |                   |Dong|BslA|Craq|Guér|H Fo|Mésa|P dC|
% +-------------------+----+----+----+----+----+----+----+                                                              
% |Donges             |   0| 148|  55|  32|  70| 140|  73|
% |Brain-s.-l’Authion | 148|   0|  93| 180|  99|  12|  72|
% |Craquefou          |  55|  93|   0|  85|  20|  83|  28|
% |Guérande           |  32|  80|  85|   0| 100| 174|  99|
% |Haie Fouassière    |  70|  99|  20| 100|   0|  85|  49|
% |Mésanger           | 140|  12|  83| 174|  85|   0|  73|
% |Ponts-de-Cé        |  73|  72|  28|  99|  49|  73|   0|
% +-------------------+----+----+----+----+----+----+----+ 
%                                                             
% The transport company uses tankers with a capacity of 39000 liters
% for the deliveries. Determine the tours for delivering to all
% clients that minimize the total number of kilometers driven.
%
% VARIABLES
%
% distance                   The distance matrix
% demand                     Demand
% capacity                   Capacity of tanker
%
% RESULTS
%
% Run this for an interpretation of the results
% Result = heatingoildeliveryEx(2);
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

function Result = heatingoildeliveryEx(PriLev)

if nargin < 1
   PriLev = 1;
end

distance     = [   0 148  55  32  70 140  73;...
      148   0  93 180  99  12  72;...
      55  93   0  85  20  83  28;...
      32 180  85   0 100 174  99;...
      70  99  20 100   0  85  49;...
      140  12  83 174  85   0  73;...
      73  72  28  99  49  73   0];

demand       = [14000;3000;6000;16000;15000;5000];
capacity     = 39000;

Prob = heatingoildelivery(distance, demand, capacity);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   s      = 7;                                         % number of sites
   names  = ['Dong'; 'BslA'; 'Craq'; 'Guér'; 'H Fo'; 'Mésa'; 'P dC'];
   tour   = reshape(Result.x_k(1:s*s),s,s);            % extract tour
   tour(find(tour<0.5)) = 0;                           % remove false zeros
   first  = find(tour(1,:));                           % might be an array
   disp(['min distance of ' num2str(Result.f_k) ' km by using'])
   
   for init = 1:length(first),                         % there might be many first stops
      this_tour = [names(1,:)];                        % start in Dong   
      site = first(init);                              
      this_tour = [this_tour ' -> ' names(site,:)];    % add city 2
      loop_me = 1;
      next = find(tour(site,:));                       % what city is next?
      
      if next == 1,                                    % if we are going home:
         loop_me = 0;                                  % do not enter while-loop
         this_tour = [this_tour  ' -> ' names(1,:)];   % just add home and quit
         disp(['Tour ' num2str(init) ': ' num2str(this_tour) ])
      end                                                     
      
      while loop_me == 1,                              % if more than one stop
         this_tour = [this_tour  ' -> ' names(next,:)];% add them one at a time
         next = find(tour(next,:));
         
         if next == 1,                                  % when we are going home
            loop_me = 0;                                % stop loop
            this_tour = [this_tour  ' -> ' names(1,:)]; % finish names and quit
            disp(['  tour ' num2str(init) ': ' num2str(this_tour) ])
         end
      end
   end
end

% MODIFICATION LOG
%
% 051020 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end

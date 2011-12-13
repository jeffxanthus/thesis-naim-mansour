%% Heating Oil Delivery
%
%% Problem description
% A transporter has to deliver heating oil from the refinery at
% Donges to a certain number of clients in the west of France. His
% clients are located at Brain-sur-l’Authion, Craquefou, Guérande, la
% Haie Fouassière, Mésanger and les Ponts-de-Cé. The following table
% lists the demands in liters for the different sites.
%
% Demands by clients (in liters)
%
%  +-----------+-----+-----+-----+-----+
%  | BslA| Craq| Guér| H Fo| Mésa| P dC|
%  +-----+-----+-----+-----+-----+-----+
%  |14000| 3000| 6000|16000|15000| 5000|
%  +-----+-----+-----+-----+-----+-----+
%
% The next table contains the distance matrix between the clients and
% the refinery.
%
% Distance matrix (in km)
%
%  +-------------------+----+----+----+----+----+----+----+
%  |                   |Dong|BslA|Craq|Guér|H Fo|Mésa|P dC|
%  +-------------------+----+----+----+----+----+----+----+
%  |Donges             |   0| 148|  55|  32|  70| 140|  73|
%  |Brain-s.-l’Authion | 148|   0|  93| 180|  99|  12|  72|
%  |Craquefou          |  55|  93|   0|  85|  20|  83|  28|
%  |Guérande           |  32|  80|  85|   0| 100| 174|  99|
%  |Haie Fouassière    |  70|  99|  20| 100|   0|  85|  49|
%  |Mésanger           | 140|  12|  83| 174|  85|   0|  73|
%  |Ponts-de-Cé        |  73|  72|  28|  99|  49|  73|   0|
%  +-------------------+----+----+----+----+----+----+----+
%
% The transport company uses tankers with a capacity of 39000 liters
% for the deliveries. Determine the tours for delivering to all
% clients that minimize the total number of kilometers driven.
%
%% Variables
%
%  distances                   The distance matrix
%  demand                     Demand
%  capacity                   Capacity of tanker
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
distances     = [   0 148  55  32  70 140  73;...
    148   0  93 180  99  12  72;...
    55  93   0  85  20  83  28;...
    32 180  85   0 100 174  99;...
    70  99  20 100   0  85  49;...
    140  12  83 174  85   0  73;...
    73  72  28  99  49  73   0];

demand       = [14000;3000;6000;16000;15000;5000];
capacity     = 39000;

n = size(distances,1);
prec = tom('prec',n,n,'int');
q = tom('q',n-1,1);

% Bounds
bnds = {0 <= prec <= 1, prec(1:n+1:n^2) == 0, q >= 0, ...
    demand <= q <= capacity*ones(n-1,1)};

% Customer constraint2.
con1 = {sum(prec(:,2:end),1) == 1};
con2 = {sum(prec(2:end,:),2) == 1};

% Quantity, demand, capacity constr.
con3 = {q <= capacity+(demand-capacity).*prec(1,2:end)'};

% Quantity, demand, capacity constr.
con4 = [];
count = 1;
for i=2:n
    for j=2:n
        if i~=j
            con4{count} = {q(j-1) >= q(i-1) + demand(j-1) - ...
                capacity + capacity*prec(i,j) + ...
                (capacity - demand(j-1) - demand(i-1))*prec(j,i)};
            count = count + 1;
        end
    end
end

% Objective
objective = sum(sum(distances.*prec));

constraints = {bnds, con1, con2, con3, con4};
options = struct;
options.solver = 'cplex';
options.name   = 'Heating Oil Delivery';
sol = ezsolve(objective,constraints,[],options);

f_k = subs(objective,sol);

PriLev = 1;
if PriLev > 0
    s      = 7;                      % number of sites
    names  = ['Dong'; 'BslA'; 'Craq'; 'Guér'; 'H Fo'; 'Mésa'; 'P dC'];
    tour   = sol.prec;               % extract tour
    tour(find(tour<0.5)) = 0;        % remove false zeros
    first  = find(tour(1,:));        % might be an array
    disp(['min distance of ' num2str(f_k) ' km by using'])

    for init = 1:length(first),      % there might be many first stops
        this_tour = [names(1,:)];    % start in Dong
        site = first(init);
        % add city 2
        this_tour = [this_tour ' -> ' names(site,:)];
        loop_me = 1;
        next = find(tour(site,:));   % what city is next?

        if next == 1,                % if we are going home:
            loop_me = 0;             % do not enter while-loop
            % just add home and quit
            this_tour = [this_tour  ' -> ' names(1,:)];
            disp(['Tour ' num2str(init) ': ' num2str(this_tour) ])
        end

        while loop_me == 1,          % if more than one stop
            % add them one at a time
            this_tour = [this_tour  ' -> ' names(next,:)];
            next = find(tour(next,:));

            if next == 1,            % when we are going home
                loop_me = 0;         % stop loop
                % finish names and quit
                this_tour = [this_tour  ' -> ' names(1,:)];
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
% 090316 med   Converted to tomSym
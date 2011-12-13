%% Planning a Flight Tour
%
%% Problem description
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
%  +--+---+---+---+---+---+---+
%  |  | A2| A3| A4| A5| A6| A7|
%  +--+---+---+---+---+---+---+
%  |A1|786|549|657|331|559|250|
%  |A2|668|979|593|224|905|   |
%  |A3|316|607|472|467|   |   |
%  |A4|890|769|400|   |   |   |
%  |A5|386|559|   |   |   |   |
%  |A6|681|   |   |   |   |   |
%  +--+---+---+---+---+---+---+
%
%% Variables
%
%  distances                   The distance matrix
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
distances = [   0   786   549    657    331    559   250;...
    786     0   668    979    593    224   905;...
    549   668     0    316    607    472   467;...
    657   979   316      0    890    769   400;...
    331   593   607    890      0    386   559;...
    559   224   472    769    386      0   681;...
    250   905   467    400    559    681     0];

n = size(distances,1);
fly = tom('fly',n,n,'int');

% All variables are binary.
bnds = {0 <= fly <= 1, fly(1:n+1:n^2) == 0};

% Cycle constraints
con = {sum(fly,1) == 1; sum(fly,2) == 1};

% Objective
objective = sum(sum(distances.*fly));

constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Planning a Flight Tour';

stop = 0;

while stop<100
    sol = ezsolve(objective,constraints,[],options);
    % Find a subtour
    idx = zeros(n,n);
    idx(1,1:n) = 1:7;
    for i=1:n
        for j=2:n
            idx(j,i) = find(sol.fly(:,idx(j-1,i)) ~= 0);
            if idx(j,i) == i
                break
            end
        end
    end
    shortmat = spones(idx);
    len = full(sum(shortmat,1));
    [val, idx2] = min(len);
    idx3 = idx(1:val, idx2);

    idx4 = zeros(val-1,1);
    for k=1:val-1
        idx4(k) = idx3(k)+n*(idx3(k+1)-1);
    end
    constraints = {constraints; sum(fly(idx4)) <= 1};

    if val == n
        break
    end
    stop = stop + 1;
end

PriLev = 1;
if PriLev > 0
    temp           = sol.fly;                    % reshape results
    not_home_again = true;                       % help to loop
    current_pos    = 1;                          % we start in 1
    route          = [];                         % initially empty
    route          = [route current_pos];        % add position
    while not_home_again,                        % loop if not home
        current_pos = find(temp(current_pos,:)); % go to next pos
        route = [route current_pos];             % add to route
        if current_pos == 1                      % if home
            not_home_again = false;              % stop loop
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
% 090308 med   Converted to tomSym
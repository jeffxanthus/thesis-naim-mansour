%% Location of Income Tax Offices
%
%% Problem description
% The income tax administration is planning to restructure the
% network of income tax offices in a region. The graph in the figure
% below shows the cities in the region and the major roads. The
% numbers within () close to the cities indicate the population in
% thousands of inhabitants. The arcs are labeled with the distances
% in kilometers. The income tax administration has determined that
% offices should be established in three cities to provide sufficient
% coverage. Where should these offices be located to minimize the
% average distance per inhabitant to the closest income tax office?
%
% Graph of towns and roads of the region
%
%
%    (15)   (10)   (12)   (18)
%    1 -15- 2 -22- 3 -18- 4
% 
%     | \       / |     |
%     |  24   16  |     |
%    18   \   /   |     |
%     |           20    12
%     |     5(5)  |     |
%     |           |     |
%     |     | \   |     |
%     |    12  24 |     |
%     |     |   \ |     |
% 
%    7 -15- 8 -30- 9 -12- 6 (24)
%  (11)    (16)   (13)
%     |     |   /  |   /
%     22   25  19  19 22
%     |     | /    | /
% 
%   10 -19- 11 -21- 12 (20)
%   (22)    (19)
%
%% Variables
%
%  population           Population of each town
%  numloc               Number of offices to start
%  lengths              The length of the roads
%  in/out               A road i goes between towns
%                       in(i) and out(i)
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
population = [15 10 12 18 5 24 11 16 13 22 19 20]';
numloc     = 3;
in         = [1 1 1 2 3 3 3 4 5 5 6  6 7  7 8  8  9  9 10 11]';
out        = [2 5 7 3 4 5 9 6 8 9 9 12 8 10 9 11 11 12 11 12]';
lengths    = [15 24 18 22 18 16 20 12 12 24 12 22 15 22 30 ...
    25 19 19 19 21]';

n1   = length(unique([in;out])); %Number of cities
n2   = length(in);
% Calculate distance matrix
d = inf*ones(n1,n1);
for i=1:n1
    d(i,i) = 0;
end
for i=1:n2
    d(in(i), out(i)) = lengths(i);
    d(out(i), in(i)) = lengths(i);
end
for i=1:n1 %b
    for j=1:n1 %c
        for k=1:n1 %d
            if j<k
                if d(j,k) > d(j,i)+d(i,k);
                    d(j,k) = d(j,i)+d(i,k);
                    d(k,j) = d(j,i)+d(i,k);
                end
            end
        end
    end
end

c = length(unique([in;out])); %Number of cities
dep = tom('dep',c,c,'int');
build = tom('build',c,1,'int');

% All variables are binary
bnds = {0 <= dep <= 1, 0 <= build <= 1};

% Building constraint
con1 = {sum(build) == numloc};

% Dependencies constraint
con2 = {sum(dep,2) == 1};

% Reality constraint
con3 = {dep <= repmat(build',c,1)};

% Objective
objective = sum(sum(repmat(population,1,c).*d.*dep));
constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Location of Income Tax Offices';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    cities = length(population);
    temp   = sol.build;
    build  = find(temp);
    goto   = sol.dep';
    disp(['Build the offices in towns ' num2str(build') ' and let'])
    for i = 1:length(build),
        disp(['   people from ' num2str(find(goto(build(i),:))) ...
            ' travel to ' num2str(build(i)) ])
    end
end

% MODIFICATION LOG
%
% 051206 med   Created
% 060118 per   Added documentation
% 060125 per   Moved disp to end
% 090325 med   Converted to tomSym
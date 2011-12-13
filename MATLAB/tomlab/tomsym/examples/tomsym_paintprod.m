%% Paint Production
%
%% Problem description
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
%  +-+--+--+--+--+--+
%  | | 1| 2| 3| 4| 5|
%  +-+--+--+--+--+--+
%  |1| 0|11| 7|13|11|
%  |2| 5| 0|13|15|15|
%  |3|13|15| 0|23|11|
%  |4| 9|13| 5| 0| 3|
%  |5| 3| 7| 7| 7| 0|
%  +-+--+--+--+--+--+
%
% Since the company also has other activities, it wishes to deal
% with this weekly production in the shortest possible time (blending
% and cleaning). Which is the corresponding order of paint batches?
% The order will be applied every week, so the cleaning time between
% the last batch of one week and the first of the following week
% needs to be counted for the total duration of cleaning.
%
%% Variables
%
%  cleantimes        Times to clean from batch i to j
%  prodtimes         Production times per batch
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
cleantimes = [ 0 11  7 13 11;...
    5  0 13 15 15;...
    13 15  0 23 11;...
    9 13  5  0  3;...
    3  7  7  7  0];

prodtimes  = [40;35;45;32;50];

n = size(cleantimes,1);
succ = tom('succ',n,n,'int');
y = tom('y',n,1);

% All slots are integers
bnds = {0 <= succ <= 1, y >= 0, succ((1:n+1:n^2)) == 0};

% Only one transition at a given time
con1 = {sum(succ,1) == 1};

% Only one transition to a given batch
con2 = {sum(succ,2) == 1};

% Sub-cycle constraint
con3 = cell(n*(n-1),1);
for i=1:n
    for j=2:n
        con3{(i-1)*(n-1)+j-1} = {y(j) >= y(i)+1-n*(1-succ(i,j))};
    end
end

% Objective
objective = sum(sum((repmat(prodtimes,1,n) + cleantimes).*succ));

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Paint Production';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    temp1 = [sol.succ, sol.y];
    link    = [];  % connections

    for i = 1:n
        for j = 1:n
            if temp1(i,j) > 0.999
                link = [[i j ]; link ]; % finding connections
            end
        end
    end

    first = link(1:1);    % start batch
    next =  link(1,2);    % next batch
    order = first;        % ordered batches

    for k = 1:n
        order = [order next];     % adding next
        next =  link(find(link(:,1)==next),2); % finding new next
    end
    disp(['one best order: ' num2str(order)])  % display solution
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060111 per   Added documentation.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
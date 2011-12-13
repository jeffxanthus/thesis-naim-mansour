%% Combining different modes of transport
%
%% Problem description
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
%   Pairs of cities
%  +----+---+---+---+---+
%  |    |1–2|2–3|3–4|4–5|
%  +----+---+---+---+---+
%  |Rail| 30| 25| 40| 60|
%  |Road| 25| 40| 45| 50|
%  |Air | 40| 20| 50| 45|
%  +----+---+---+---+---+
%
% The next table summarizes the costs for changing the mode of
% transport in $ per tonne. The cost is independent of location.
%
% Cost for changing the mode of transport
%
%  +---------+----+----+---+
%  |from \ to|Rail|Road|Air|
%  +---------+----+----+---+
%  |Rail     |  0 |  5 | 12|
%  |Road     |  8 |  0 | 10|
%  |Air      | 15 | 10 |  0|
%  +---------+----+----+---+
%
% How should we organize the transport of the load at the least cost?
%
%% Variables
%
%  transpcost                 Transport costs
%  changecost                 Cost to change mode of transport
%  demand                     Load to transport
%
%% References
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 11, 2009.

%% Problem setup
transpcost = tomArray([ 30 25 40 60 ;...
    25 40 45 50 ;...
    40 20 50 45]);

changecost = tomArray([  0  5  12;...
    8  0  10;...
    15 10   0]);

demand = 30;

m = size(changecost,1); % mode
l = size(transpcost,2); % legs
n = size(changecost,2); % mode

i  = tomArrayIdx('i',1:m);    % index over all modes to change from
j  = tomArrayIdx('j',1:l);    % index over all legs
jc = tomArrayIdx('jc',1:l-1); % index over all legs that involve a change
k  = tomArrayIdx('k',1:n);    % index over all modes to chang to

change = tom('change',m*n*(l-1),1,'int');
change = tomArray(change,[m,n,l-1]);

use = tom('use', m*l, 1, 'int');
use = tomArray(use,[m,l]);

% All variables are binary.
bnds = {0 <= change <= 1, 0 <= use <= 1};

% Single mode of transport
con1 = {sum(use(i,j),i) == 1};

% Single change of mode
con2 = {sum(sum(change(i,k,jc),i),k) == 1};

% Change to mode relationship
con3 = use(i,jc) + use(k,jc+1) >= 2*change(i,k,jc);

% Objective
objective = sum(vec(transpcost(i,j).*use(i,j))) + ...
    sum(vec(changecost(i,k).*change(i,k,jc)));

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Combining Diff Modes of Transp';
[sol, Result] = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    tmodes  = reshape(sol.use,3,4);
    means  = ['rail';'road';'air '];
    disp(['the min cost (' num2str(Result.f_k) ') is found by,'])
    for m = 1:size(tmodes,2),
        disp(['   going from town ' num2str(m) ' to town ' num2str(m+1) ...
            ' by ' means(find(tmodes(:,m)),:)])
    end
end

% MODIFICATION LOG
%
% 051020 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end
% 090407 med   Converted to tomSym
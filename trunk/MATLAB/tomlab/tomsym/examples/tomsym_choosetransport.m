%% Choosing the Mode of Transport
%
%% Problem description
% A company in the south-west of France needs to transport 180 tonnes
% of chemical products stored in depots D1 to D4 to the three
% recycling centers C1, C2, and C3. The depots D1 to D4 contain
% respectively 50, 40, 35, and 65 tonnes, that is 190 tonnes in
% total. Two modes of transport are available: road and rail. Depot
% D1 only delivers to centers C1 and C2 and that by road at a cost of
% $ 12k/t and $ 11k/t. Depot D2 only delivers to C2, by rail or road
% at $ 12k/t and $ 14k/t respectively. Depot D3 delivers to center C2
% by road ($ 9k/t) and to C3 by rail or road for $ 4k/t and $ 5k/t
% respectively. The depot D4 delivers to center C2 by rail or road at
% a cost of $ 11k/t and $ 14k/t, and to C3 by rail or road at $ 10k/t
% and $ 14k/t respectively.
%
% Its contract with the railway company for the transport of chemical
% products requires the company to transport at least 10 tonnes and
% at most 50 tonnes for any single delivery. Besides the standard
% security regulations, there are no specific limitations that apply
% to road transport. How should the company transport the 180 tonnes
% of chemicals to minimize the total cost of transport?
%
% Graph of connections.
%
%                    ------  7 (road) \
%       +--  2 (D1)        /           \
%      /            \     /             12 (C1)
%     /              \   /             /
%    /              --+--    6 (rail) /         \
%   |    +-  3 (D2)   |     /                    \
%   |   /           --+----/                      \
%   |  /              |                            \
%   | /               |                             \
%                +----+----  8 (rail) \              \
%   1           /      \               \
%     \        /        \               13 (C2) ----- 15
%   |  \      /          \             /
%   |   \                 \  9 (road) /              /
%   |    +-- 5 ------------                         /
%   |     (D4) \     ------                        /
%   |           \   /                             /
%    \           +--+-----  10 (rail) \          /
%     \             /     /            \        /
%      \       ----+     /              14 (C3)
%       +--- 4 ---------+              /
%        (D3) \             11 (road) /
%              +-----------
%
%
%% Variables
%
%  arcs_out/in                For an arc i arcs_out(i) is its starting node
%                             and arcs_in(i) ts target node
%  arcs_min                   Minimal load for an arc
%  arcs_max                   Maximal load for an arc
%  arcs_cost                  Cost for an arc
%  minflow                    Flow from 1 to 15
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Oct 7, 2005.   Last modified Mar 8, 2009.

%% Problem setup
arcs_out        = [ 1  1  1  1  2  2  3  3  4  4  4  5  5  5  5 ...
    6  7  8  9 10 11 12 13 14 ]';
arcs_in         = [ 2  3  4  5  7  9  6  7  9 10 11  8  9 10 11 ...
    12 12 13 13 14 14 15 15 15 ]';
arcs_min        = [zeros(1,6) 10 0 0 10 0 10 0 10 0 ...
    zeros(1,9)]';
arcs_max        = [50 30 35 65 inf inf 50 inf inf 50 inf 50 inf 50 inf ...
    inf*ones(1,9)]';
arcs_cost       = [0 0 0 0 12 11 12 14 9 4 5 11 14 10 14 ...
    zeros(1,9) ]';
minflow         = 50+30+35+65;

n  = length(arcs_out); % arcs
n1 = length(unique([arcs_in ; arcs_out]));

arcs = tom('arcs',n,1,'int');

% Bounds
bnds = {arcs_min <= arcs <= arcs_max};

% Node constraint, except for source and sink
con1 = cell(n1-2,1);
for i=2:n1-1
    idx_in  = find(arcs_in  == i);
    idx_out = find(arcs_out == i);
    con1{i-1} = {sum(arcs(idx_in)) == sum(arcs(idx_out))};
end

% Source constraint
con2 = {sum(arcs(arcs_out == 1)) == minflow};

% Objective
objective = arcs_cost'*arcs;

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Choosing the Mode of Transport';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    arcs = sol.arcs;
    disp('transport as follows:')
    for i = 1:length(arcs),
        if arcs(i) ~= 0,
            disp(['    ' num2str(arcs(i)) ' units from node ' ...
                num2str(arcs_out(i)) ' to node ' num2str(arcs_in(i))])
        end
    end
end

% MODIFICATION LOG
%
% 051019 med   Created
% 060112 per   Added documentation
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
%% Routing Telephone Calls
%
%% Problem description
% A private telephone company exploits a network, represented in the
% figure below, between five cities: Paris, Nantes, Nice, Troyes, and
% Valenciennes.
%
% Structure of the network of the company
%
%  Nantes --(300)-- Paris --(200)-- Valenciennes
% 
%      \            /                |
%     (120)      (300)              (70)
%        \        /                  |
% 
%           Nice   ----(80)----  Troyes
%
% The number beside each edge (connection) is the capacity of the
% link in terms of circuits. This issue is worth some further
% explanation. Suppose that a person A at Nantes calls a person B at
% Nice. The company needs to find a path formed by non-saturated
% edges from Nantes to Nice to route the binary flow corresponding to
% the digitized voice of A. But the conversation is only possible if
% a flow from Nice to Nantes is established so that B may answer A.
% In digital telephone networks, the binary flows all have the same
% throughput standard (often 64 kbps). The associated capacity is
% called a channel. The return channel uses the same edges as the
% first channel, but in the opposite direction. This linked pair of
% channels necessary for a telephone conversation is a circuit.
%
% The routing of the return channel through the same edges is due to
% the fact that the call could fail if one waited until the callee
% picked up the phone before searching for a non-saturated return
% path. This is why, at the moment of the call, the routing system
% constructs the path edge by edge and immediately reserves the edges
% of the return channel. As a consequence, as the capacity of an edge
% is consumed in increments of a bidirectional circuit, we do not
% consider any directed flows. For example, there may be 10 persons
% calling from Nantes to Nice and 20 from Nice to Nantes, or the
% opposite: in both cases, 30 circuits are used.
%
% At a given moment, the company is facing demands for circuits given
% in the following table. Is it possible to satisfy the demands
% entirely? If this is not possible, try to transmit as much as
% possible. In every case, indicate the corresponding routing,
% that is, the access paths used.
%
% Demand of circuits
%
%  +-------------------+--------+
%  |Cities             |Circuits|
%  +-------------------+--------+
%  |Nantes-Nice        |  100   |
%  |Nantes-Troyes      |   80   |
%  |Nantes-Valenciennes|   75   |
%  |Nice-Valenciennes  |  100   |
%  |Paris-Troyes       |   70   |
%  +-------------------+--------+
%
%% Variables
%
%  arcs_out/in          For an arc i arcs_out(i) is its starting node
%                       and arcs_in(i) ts target node
%  capacity             Capacity of arcs
%  demand_out/in        For a demand i demand_out(i) is its starting node
%                       and demand_in(i) its target node
%  demands              Demand
%  path_mat             Possible paths
%  paths                In/Out of possible paths
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
% Nantes, Paris, Nice, Valenciennes, Troyes

arcs_in    = [1 1 2 2 3 4]';
arcs_out   = [2 3 3 4 5 5]';
capacity   = [300 120 300 200 80 70]';

demand_in  = [1 1 1 3 2]';
demand_out = [3 5 4 4 5]';
demands    = [100 80 75 100 70]';

path_mat    = [1 3 0 0 0;...
    1 2 3 0 0;...
    1 2 4 5 3;...
    1 2 4 5 0;...
    1 2 3 5 0;...
    1 3 5 0 0;...
    1 3 2 4 5;...
    1 2 4 0 0;...
    1 3 2 4 0;...
    1 2 3 5 4;...
    1 3 5 4 0;...
    3 1 2 4 0;...
    3 2 4 0 0;...
    3 5 4 0 0;...
    2 4 5 0 0;...
    2 1 3 5 0;...
    2 3 5 0 0];

paths = [1 3;...
    1 3;...
    1 3;...
    1 5;...
    1 5;...
    1 5;...
    1 5;...
    1 4;...
    1 4;...
    1 4;...
    1 4;...
    3 4;...
    3 4;...
    3 4;...
    2 5;...
    2 5;...
    2 5];

n  = size(paths,1);   % paths

flows = tom('flows',n,1);

% No variables are binary.
bnds = {0 <= flows};

% Capacity constraint on arc
n1 = length(arcs_in);
con1 = cell(n1,1);
for i=1:n1
    [idx_i, idx_j] = find(path_mat(1:end,1:end-1) == arcs_in(i));
    next_nodes = path_mat(idx_i+n*idx_j);
    idx2 = find(next_nodes == arcs_out(i));
    new_idx = idx_i(idx2);
    con1{i} = {sum(flows(new_idx)) <= capacity(i)};
end

% Demand constraint
n2 = length(demands);
con2 = cell(n2,1);
for i=1:n2
    [idx_i, idx_j] = find(paths(1:end,1) == demand_in(i));
    next_nodes = paths(idx_i+n*idx_j);
    idx2 = find(next_nodes == demand_out(i));
    new_idx = idx_i(idx2);
    con2{i} = {sum(flows(new_idx)) <= demands(i)};
end

% Objective
objective = -sum(flows);

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Routing Telephone Calls';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    paths  = sol.flows;
    idx = find(paths);
    load = [];
    cities = ['Nant'; 'Pari'; 'Nice'; 'Vale'; 'Troy'];
    disp('INCOMING AND OUTGOING CIRCUITS')
    for i = 1:length(idx),
        id = idx(i);
        temp_path = path_mat(id,:);
        temp_path = temp_path(find(temp_path));
        city_path = [];
        for city = 1:length(temp_path),
            city = temp_path(city);
            city_path = [city_path cities(city,:) ' ' ];
        end
        disp(['   ' num2str(paths(id)) ' circuits ' ...
            cities(temp_path(1),:) '-' cities(temp_path(end),:) ...
            ' using the following path: ' num2str(city_path)])
        for j = 2:length(temp_path),
            out = temp_path(j-1);
            in  = temp_path(j);
            if in < out,
                in_temp = in;
                in = out;
                out = in_temp;
            end
            if (size(load) > [out in]),
                load(out,in) = load(out,in) + paths(id);
            else
                load(out,in) = paths(id);
            end
        end
    end
    disp('LOAD BETWEEN CITIES')
    for i = 1:length(cities)-1,
        for j = 1:length(cities),
            if load(i,j) > 0,
                disp(['   between ' cities(i,:) ' and ' cities(j,:) ': ' ...
                    num2str(load(i,j)) ' circuits.' ])
            end
        end
    end
end

% MODIFICATION LOG
%
% 051122 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
%% Water Supply Management
%
%% Problem description
% The graph in the figure below shows a water transport network. The
% nodes, numbered from 1 to 10, represent the cities, the reservoirs,
% and the pumping stations connected by a network of pipes. The three
% cities Gotham City, Metropolis, and Spider Ville are supplied from
% two reservoirs. The availabilities of water from these reservoirs
% in thousands of m3/h are 35 for reservoir 1 and 25 for reservoir 2.
% The capacity of each pipe connection is indicated in thousands of
% m3/h in the graph.
%
% A study is undertaken to find out whether the existing network will
% be able to satisfy the demands of the cities in ten years time,
% that is 18, 15, and 20 thousand m3/h. Determine the maximum flow in
% the current network. Will it be sufficient in ten years from now?
%
% Water transport network
%
%
%     35          20        15         7          18
%    Reservoir-1 ----> 3 -------> 4 ------->  8-Gotham city
%                                              ^
%              | \15   |            \       ---+
%              |  \    |10           ----\ /
%            12|   \   |                  X
%               \   \  |       10        / \10
%                \   V V  --------------/   ---+
%     25          \      /      15             V   15
%    Reservoir-2 --+-> 5 ------------------->  9-Metropolis
%                  | 6   \       15            ^
%                \ \      -----------   -------+
%                 \ \                \ /   10
%                  \ \12              X
%                22 \ \              / \
%                    V V            /   -----+
%                          22         10     V     20
%                      6 -------> 7 -------> 10-Spider Ville
%
%
%% Variables
%
%  arcs_out/in       For an arc i arcs_out(i) is its target node
%                    and arcs_in(i) its starting node
%  capacity          Capacity of an arc
%  source            Artificial node acting as a source
%  sink              Artificial node acting as a sink
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

arcs_in  = [ 1  1  1  2  2  3  3  4  4  5  5  5  6  7  7  8  9 10 11 11]';
arcs_out = [ 3  5  6  5  6  4  5  8  9  8  9 10  7  9 10 12 12 12  1  2]';

capacity = [20 15 12  6 22 15 10  7 10 10 15 15 22 10 10 18 15 20 35 25]';
source   = 11;
sink     = 12;

n = length(arcs_in);  %arcs

arcs = tom('arcs',n,1,'int');

% All variables are binary
bnds = {0 <= arcs <= capacity};

% Kirchhoffs's law
n1   = length(unique([arcs_in;arcs_out])) - 2;
b_L  = zeros(n1,1);
b_U  = zeros(n1,1);
con = cell(n1,1);
count = 1;
for i=1:n1+2
    if ~(i == source || i == sink)
        idx1 = find(arcs_in == i);
        idx2 = find(arcs_out == i);
        con{count} = {sum(arcs(idx1)) == sum(arcs(idx2))};
        count = count + 1;
    end
end

% Objective
idx = find(arcs_out == sink);
objective = -sum(arcs(idx));
constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Water Conveyance';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    x        = sol.arcs;
    names    = ['Reservoir 1 '; ...
        'Reservoir 2 '; ...
        'Node 3      '; ...
        'Node 4      '; ...
        'Node 5      '; ...
        'Node 6      '; ...
        'Node 7      '; ...
        'Gotham City '; ...
        'Metropolis  '; ...
        'Spider Ville'; ...
        'Production  '; ... % this is the source node
        'Consumption '];    % this is the sink node
    disp('Maximum flow of the network is as follows:')
    for i = 1:length(arcs_in),
        if x(i) ~= 0,
            disp([names(arcs_in(i),:) ' -> ' names(arcs_out(i),:) ': ' num2str(x(i)) ])
        end
    end
end

% MODIFICATION LOG
%
% 051205 med   Created.
% 060117 per   Added documentation.
% 060125 per   Moved disp to end
% 090325 med   Converted to tomSym
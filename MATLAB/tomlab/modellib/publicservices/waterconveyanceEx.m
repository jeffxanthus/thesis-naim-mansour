% function Result = waterconveyanceEx(PriLev)
%
% Creates a TOMLAB MILP problem for water conveyance
%
% WATER CONVEYANCE / WATER SUPPLY MANAGEMENT
%
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
% VARIABLES
%
% arcs_out/in                For an arc i arcs_out(i) is its target node
%                            and arcs_in(i) its starting node
% capacity                   Capacity of an arc
% source                     Artificial node acting as a source
% sink                       Artificial node acting as a sink
%
% RESULTS
%
% For an interpretation of the results, run:
% Result   = waterconveyanceEx(2);
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
% Result       Result structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Result = waterconveyanceEx(PriLev)

if nargin < 1
   PriLev = 1;
end

arcs_in         = [ 1  1  1  2  2  3  3  4  4  5  5  5  6  7  7  8  9 10 11 11]';
arcs_out        = [ 3  5  6  5  6  4  5  8  9  8  9 10  7  9 10 12 12 12  1  2]';

capacity        = [20 15 12  6 22 15 10  7 10 10 15 15 22 10 10 18 15 20 35 25]';
source          = 11;
sink            = 12;

Prob = waterconveyance(arcs_in, arcs_out, capacity, source, sink);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   x        = Result.x_k;
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

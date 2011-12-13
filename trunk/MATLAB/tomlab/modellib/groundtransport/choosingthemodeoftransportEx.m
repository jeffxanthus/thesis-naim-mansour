% function Result = choosingthemodeoftransportEx(PriLev)
%
% Creates a TOMLAB MIP problem for choosing the mode of transport
%
% CHOOSING THE MODE OF TRANSPORT
%
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
% VARIABLES
%
% arcs_out/in                For an arc i arcs_out(i) is its starting node
%                            and arcs_in(i) ts target node
% arcs_min                   Minimal load for an arc
% arcs_max                   Maximal load for an arc
% arcs_cost                  Cost for an arc
% minflow                    Flow from 1 to 15
%
% RESULTS
%
% For an interpretation of the results, try the following
% Result   = choosingthemodeoftransportEx(2);
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
% Written Oct 19, 2005.   Last modified Oct 19, 2005.

function Result = choosingthemodeoftransportEx(PriLev)

if nargin < 1
   PriLev = 1;
end

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

Prob = choosingthemodeoftransport(arcs_out, arcs_in, arcs_min, arcs_max, arcs_cost, minflow);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   arcs     = Result.x_k;
   
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
% 051019 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end

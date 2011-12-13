% function Result = networkreliabilityEx(PriLev)
%
% Creates a TOMLAB MIP problem for network reliability
%
% NETWORK RELIABILITY
%
% We consider the military telecommunications network represented
% below. This network consists of eleven sites connected by
% bidirectional lines for data transmission. For reliability reasons
% in the case of a conflict, the specifications require that the two
% sites (nodes) 10 and 11 of the network remain able to communicate
% even if any three other sites of the network are destroyed. Does
% the network satisfy this requirement?
%
%
%         1 --- 2 ---------- 8
%
%      /  |   / |          / |
%     /   |  /  |         /  |
%    /    | /   |        /   |
%               |            |
%  11 --- 3  ---+---- 10     |
%               |            |
%   | \ /   \   |   / |  \   |
%   |  X     \  |  /  |   |  |
%   | / \     \ | /   |   |  |
%                     |   |  |
%   4 --- 5 --- 9     |   |  |
%                     |   |  |
%     \         |    /    |  |
%      \        |   /     \  |
%       \       |  /       \ |
%        \                   
%         ----- 6 ---------- 7
%
%
%
% VARIABLES
%
% arcs_out/in                For an arc i arcs_out(i) is its starting node
%                            and arcs_in(i) ts target node
%
% RESULTS
%
% for an interpretation of the results, use a PriLev > 1, for example
% Result =   networkreliabilityEx(2);
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
% Written Nov 7, 2005.   Last modified Nov 7, 2005.

function Result = networkreliabilityEx(PriLev)

if nargin < 1
   PriLev = 1;
end

arcs_in    = [1 1  1 2 2 2 3 3  3  3 4 4  4 5  5 6 6  6 7  7  8  9]';
arcs_out   = [2 3 11 3 8 9 4 9 10 11 5 6 11 9 11 7 9 10 8 10 10 10]';

Prob = networkreliability(arcs_in, arcs_out);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   f      = - Result.f_k;
   x      = - Result.x_k;
   disp(['The network can manage the loss of ' num2str(f-1) ' nodes.'])
   disp(['There are ' num2str(f) ' disjunctive paths.'])
   [arcs,direction] = find(reshape(x,length(arcs_in),2));
   disp('For a maximal number of disjunctive paths, activate...')
   for arc = 1:length(arcs),
      if direction(arc) == 1,
         disp([' arc ' num2str([arcs_out(arcs(arc))]) '->' ...
               num2str([arcs_in(arcs(arc))])])
      end
      if direction(arc) == 2,
         disp([' arc ' num2str([arcs_in(arcs(arc))]) '->'  ...
               num2str([arcs_out(arcs(arc))])])
      end
   end
end

% MODIFICATION LOG
%
% 051107 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end

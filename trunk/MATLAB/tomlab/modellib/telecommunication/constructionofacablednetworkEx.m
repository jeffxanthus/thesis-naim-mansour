% function Result = constructionofacablednetworkEx(PriLev)
%
% Creates a TOMLAB MIP problem for construction of a cabled network
%
% CONSTRUCTION OF A CABLED NETWORK
%
% A university wishes to connect six terminals located in different
% buildings of its campus. The distances, in meters, between the
% different terminals are given in the table below. 
%
% Distances between the different terminals (in meters)
% 
% +----------+---+---+---+---+---+---+
% |          | T1| T2| T3| T4| T5| T6|
% +----------+---+---+---+---+---+---+
% |Terminal 1|  0|120| 92|265|149|194|
% |Terminal 2|120|  0|141|170| 93|164|
% |Terminal 3| 92|141|  0|218|103|116|
% |Terminal 4|265|170|218|  0|110|126|
% |Terminal 5|149| 93|103|110|  0| 72|
% |Terminal 6|194|164|116|126| 72|  0|
% +----------+---+---+---+---+---+---+
%
% These terminals are to be connected via underground cables. We
% suppose the cost of connecting two terminals is proportional to the
% distance between them. Determine the connections to install to
% minimize the total cost.
%
% VARIABLES
%
% distances                  a distance matrix
%
% RESULTS
%
% for an interpretation of the results, let PriLev > 1, for example:
% Result = constructionofacablednetworkEx(2); % call solver
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
% Written Nov 22, 2005.   Last modified Nov 22, 2005.

function Result = constructionofacablednetworkEx(PriLev)

if nargin < 1
   PriLev = 1;
end

distances    = [  0 120  92 265 149 194;...
      120   0 141 170  93 164;...
      92 141   0 218 103 116;...
      265 170 218   0 110 126;...
      149  93 103 110   0  72;...
      194 164 116 126  72   0];

Prob = constructionofacablednetwork(distances);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   t      = size(distances,1);                 % number of terminals
   s      = sum(1:t-1);                        % possible connections
   arcs   = [];                                % empty set of arcs
   count  = 1;                                 % counter
   for i = 1:t-1,                              % we catch true arcs
      for j = i+1:t,
         if Result.x_k(count) == 1,
            arcs  = [arcs ; i j];                
         end
         count = count + 1;
      end
   end
   for arc = 1:length(arcs),
      arc = arcs(arc,:);
      disp(['connect the terminals ' num2str(arc(1)) ' and ' num2str(arc(2)) ]) 
   end
end

% MODIFICATION LOG
%
% 051122 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end

% function Result = cuttingsteelbarsfordesklegsEx(PriLev)
%
% Creates a TOMLAB MIP problem for cutting steel bars for desk legs.
%
% CUTTING STEEL BARS FOR DESK LEGS
%
% The company SchoolDesk produces desks of different sizes for
% kindergartens, primary and secondary schools, and colleges. The
% legs of the desks all have the same diameter, with different
% lengths: 40 cm for the smallest ones, 60 cm for medium height, and
% 70 cm for the largest ones. These legs are cut from steel bars of
% 1.5 or 2 meters. The company has received an order for 108 small,
% 125 medium and 100 large desks. How should this order be produced
% if the company wishes to minimize the trim loss? 
%
% VARIABLES
%
% patterns                   The different cutting patterns
% demand                     Demand of the different lengths
% loss                       Loss for each pattern
% lengths                    The lengths
%
% RESULTS
%
% To interpret the results run this:
% Result = cuttingsteelbarsfordesklegsEx(2);
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
% Copyright (c) 2005-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Oct 10, 2005.   Last modified Dec 18, 2007.

function Result = cuttingsteelbarsfordesklegsEx(PriLev)

if nargin < 1
   PriLev = 1;
end

patterns = [0 0 2 0 2 3 0 0 1 3 0 5;...
      0 1 0 2 1 0 1 2 0 0 3 0;...
      2 1 1 0 0 0 2 1 2 1 0 0];

demand   = [108;125;100]*4;
loss     = [10 20  0 30 10 30  0 10 20 10 20  0]';
lengths  = [150;200];

Prob = cuttingsteelbarsfordesklegs(demand, patterns, lengths);
Prob.fConstant = -75280; % Constant to deduct from objective
Result = tomRun('cplex', Prob, PriLev);
Result.loss = sum(Result.x_k(:,1).*loss);

if PriLev > 1,
   x      = Result.x_k;
   idx    = find(x);
   order  = [x(idx) idx];
   disp(['a minimal loss of ' num2str(Result.loss) ' is found with this combination:' ])
   for i = 1:length(idx),
      disp([' cut ' num2str([order(i,1)]) ' bar(s) in pattern ' num2str([order(i,2)])])
   end  
end

% MODIFICATION LOG
%
% 051010 med   Created
% 051208 med   Loss added to Result
% 060112 per   Added documentation.
% 060112 per   Minor update.
% 060125 per   Moved disp to end
% 071218 ango  Multiple CPLEX solutions handled gracefully
% function Result = cuttingsheetmetalEx(PriLev)
%
% Creates a TOMLAB MIP problem for cutting smaller parts from big ones.
%
% CUTTING SHEET METAL
%
% A sheet metal workshop cuts pieces of sheet metal from large
% rectangular sheets of 48 decimeters × 96 decimeters (dm). It has
% received an order for 8 rectangular pieces of 36 dm × 50 dm, 13
% sheets of 24 dm × 36 dm, 5 sheets of 20 dm × 60 dm, and 15 sheets
% of 18 dm × 30 dm. Theses pieces of sheet metal need to be cut from
% the available large pieces. How can this order by satisfied by
% using the least number of large sheets?
%
% VARIABLES
%
% patterns                   Each column represent a possible sheet.
% demand                     Wanted rectangles
%
% RESULTS
%
% To explain the results, run:
% Result = cuttingsheetmetalEx(2);
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
% Written Oct 10, 2005.   Last modified Oct 10, 2005.

function Result = cuttingsheetmetalEx(PriLev)

if nargin < 1
   PriLev = 1;
end

% Have to manually evaluate all the possible patterns.

patterns = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
      2 1 0 2 1 0 3 2 1 0 5 4 3 2 1 0;...
      0 0 0 2 2 2 1 1 1 1 0 0 0 0 0 0;...
      0 1 3 0 1 3 0 2 3 5 0 1 3 5 6 8];
demand = [8;13;5;15];

Prob = cuttingsheetmetal(demand, patterns);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   x      = Result.x_k;
   idx    = find(x);
   order  = [x(idx) idx];
   for i = 1:length(idx),
      disp(['cut ' num2str([order(i,1)]) ' sheet(s) in pattern ' num2str([order(i,2)])])
   end
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end

%% Cutting Sheet Metal
%
%% Problem description
% A sheet metal workshop cuts pieces of sheet metal from large
% rectangular sheets of 48 decimeters × 96 decimeters (dm). It has
% received an order for 8 rectangular pieces of 36 dm × 50 dm, 13
% sheets of 24 dm × 36 dm, 5 sheets of 20 dm × 60 dm, and 15 sheets
% of 18 dm × 30 dm. Theses pieces of sheet metal need to be cut from
% the available large pieces. How can this order by satisfied by
% using the least number of large sheets?
%
%% Variables
%
%  patterns                   Each column represent a possible sheet.
%  demand                     Wanted rectangles
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

% Have to manually evaluate all the possible patterns.

%% Problem setup
patterns = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;...
    2 1 0 2 1 0 3 2 1 0 5 4 3 2 1 0;...
    0 0 0 2 2 2 1 1 1 1 0 0 0 0 0 0;...
    0 1 3 0 1 3 0 2 3 5 0 1 3 5 6 8];
demand = [8;13;5;15];

p = size(patterns,2);
use = tom('use',p,1,'int');

% Bounds
bnds = {0 <= use};

% Minimum demand must be met
con = {patterns*use >= demand};

% Objective
objective = sum(use);

constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Cutting Sheet Metal';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    x      = sol.use;
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
% 090308 med   Converted to tomSym
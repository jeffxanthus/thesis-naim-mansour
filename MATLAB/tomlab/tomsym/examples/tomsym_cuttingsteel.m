%% Cutting Steel Bars for Desk Legs
%
%% Problem description
% The company SchoolDesk produces desks of different sizes for
% kindergartens, primary and secondary schools, and colleges. The
% legs of the desks all have the same diameter, with different
% lengths: 40 cm for the smallest ones, 60 cm for medium height, and
% 70 cm for the largest ones. These legs are cut from steel bars of
% 1.5 or 2 meters. The company has received an order for 108 small,
% 125 medium and 100 large desks. How should this order be produced
% if the company wishes to minimize the trim loss?
%
%% Variables
%
%  patterns                   The different cutting patterns
%  demand                     Demand of the different lengths
%  loss                       Loss for each pattern
%  lengths                    The lengths
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
patterns = [0 0 2 0 2 3 0 0 1 3 0 5;...
    0 1 0 2 1 0 1 2 0 0 3 0;...
    2 1 1 0 0 0 2 1 2 1 0 0];

demand   = [108;125;100]*4;
loss     = [10 20  0 30 10 30  0 10 20 10 20  0]';
lengths  = [150;200];

p = size(patterns,2);
use = tom('use',p,1,'int');

% Bounds
bnds = {use >= 0};

% Minimum demand must be met
con = {patterns*use >= demand};

% Objective
objective = sum(lengths(1)*use(1:end/2) + lengths(2)*use(end/2+1:end));
% Constant to deduct from objective
objective = objective-75280;
constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Cutting Steel Bars';
sol = ezsolve(objective,constraints,[],options);

loss = sum(sol.use.*loss);

PriLev = 1;
if PriLev > 0
    x      = sol.use;
    idx    = find(x);
    order  = [x(idx) idx];
    disp(['a minimal loss of ' num2str(loss) ' is found with this combination:' ])
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
% 090308 med   Converted to tomSym
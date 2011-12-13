%% Planning the Production of Bicycles
%
%% Problem description
% A company produces bicycles for children. The sales forecast in
% thousand of units for the coming year are given in the following
% table. The company has a capacity of 30,000 bicycles per month.
% It is possible to augment the production by up to 50% through
% overtime working, but this increases the production cost for a
% bicycle from the usual $ 32 to $ 40.
%
% Sales forecasts for the coming year in thousand units
%
%  +---+---+---+---+---+---+---+---+---+---+---+---+
%  |Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|
%  +---+---+---+---+---+---+---+---+---+---+---+---+
%  | 30| 15| 15| 25| 33| 40| 45| 45| 26| 14| 25| 30|
%  +---+---+---+---+---+---+---+---+---+---+---+---+
%
% Currently there are 2,000 bicycles in stock. The storage costs have
% been calculated as $ 5 per unit held in stock at the end of a month.
% We assume that the storage capacity at the company is virtually
% unlimited (in practice this means that the real capacity, that is
% quite obviously limited, does not impose any limits in our case).
% We are at the first of January. Which quantities need to be
% produced and stored in the course of the next twelve months in order
% to satisfy the forecast demand and minimize the total cost?
%
%% Variables
%
%  normcapacity          the normal production capacity
%  extracapacity         extra capacity
%  normcost              normal cost
%  extracost             cost per bike if overtime
%  demand                bikes wanted per month
%  startstock            bikes in store
%  storagecost           cost to have a bike in store one month
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
normcapacity  = 30000;
extracapacity = 15000;
normcost      = 32;
extracost     = 40;
demand        = [30;15;15;25;33;40;45;45;26;14;25;30]*1000;
startstock    = 2000;
storagecost   = 5;

t = length(demand); % Months

pnorm = tom('pnorm',t,1,'int');
pover = tom('pover',t,1,'int');
store = tom('store',t,1,'int');

% All slots are integers
bnds = {0 <= pnorm <= normcapacity, 0 <= pover <= extracapacity,...
    0 <= store};

% Initial demand constraint
con1 = {pnorm(1) + pover(1) + startstock == demand(1) + store(1)};

% Monthly demand constraint
con2 = {pnorm(2:end) + pover(2:end) + store(1:end-1) == ...
    demand(2:end) + store(2:end)};

% Objective
objective = sum(normcost*pnorm) + sum(extracost*pover) + ...
    sum(storagecost*store);

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Planning the Production of Bicycles';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    Nmonths = length(demand);
    temp = [sol.pnorm, sol.pover, sol.store];
    for i = 1:Nmonths,
        disp(['Month ' num2str(i) ':'])
        disp(['   produce ' num2str(temp(i,1)) ' regular bikes' ])
        if temp(i,2) > 0,
            disp(['   and     ' num2str(temp(i,2)) ' extras' ])
        end
        if temp(i,3) > 0,
            disp(['   let     ' num2str(temp(i,3)) ' be stored' ])
        end
    end
end

% MODIFICATION LOG
%
% 051017 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
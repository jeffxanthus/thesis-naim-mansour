%% Fleet Planning for Vans
%
%% Problem description
% A chain of department stores uses a fleet of vans rented from
% different rental agencies. For the next six months period it has
% forecast the following needs for vans (table below):
%
% Requirements for vans for six months
%
%  +---+---+---+---+---+---+
%  |Jan|Feb|Mar|Apr|May|Jun|
%  +---+---+---+---+---+---+
%  |430|410|440|390|425|450|
%  +---+---+---+---+---+---+
%
% At the 1st January, the chain has 200 vans, for which the rental
% period terminates at the end of February.
%
% To satisfy its needs, the chain has a choice among three types of
% contracts that may start the first day of every month: 3-months
% contracts for a total cost of $1700 per van, 4-months contracts at
% $2200 per van, and 5-months contracts at $2600 per van. How many
% contracts of the different types need to be started every month in
% order to satisfy the company�s needs at the least cost and to have
% no remaining vans rented after the end of June?
%
%  Running contracts in month 5 (May)   .........
%                                       :       :
%                               +-------:-------:-------+
%                               | Rent34:       :       |
%                       +-------+-------:-------:-------+
%                       |     Rent33    :       :
%                       +---------------:-------:-------+
%                       | Rent43        :       :       |
%               +-------+---------------:-------:-------+
%               |    Rent42             :       :
%               +-----------------------:-------:-------+
%               | Rent52                :       :       |
%       +-------+-----------------------:-------:-------+
%       |    Rent51                     :       :
%       +-------------------------------:-------:
%               .       .       .       :       :       .
%  Month   Jan  :  Feb  :  Mar  :  Apr  :  May  :  Jun  :
%       :.......:.......:.......:.......:       :.......:
%                                       :.......:
%
%
%% Variables
%
%  demand                     Demand of vans per month
%  initialsupply              Vans per month
%  contractlength             Contracts available
%  contractcost               Cost of contracts
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
demand         = [ 430; 410; 440; 390; 425; 450];
initialsupply  = [ 200; 200;   0;   0;   0;   0];
contractlength = [ 3; 4; 5];
contractcost   = [ 1700; 2200; 2600];

c = length(contractlength); % contracts
m = length(demand);         % months

rent = tom('rent',c,m,'int');

% All variables are integers.
bnds = {0 <= rent};

% Demand constraint
for i=1:m
    for j=1:c
        idx1(i,j) = max(1,i-contractlength(j)+1);
        idx2(i,j) = min(i,m-contractlength(j)+1);
    end
end

con1 = cell(m,1);
for i=1:m
    con1{i} = 0;
    for j=1:c
        con1{i} = con1{i} + sum(rent(j,idx1(i,j):idx2(i,j)));
    end
    con1{i} = {con1{i} >= demand(i) - initialsupply(i)};
end

% Objective
objective = sum(contractcost'*rent);

constraints = {bnds, con1};
options = struct;
options.solver = 'cplex';
options.name   = 'Fleet Planning for Vans';
sol = ezsolve(objective,constraints,[],options);

f_k = subs(objective,sol);

PriLev = 1;
if PriLev > 0
    Nmonths = length(demand);
    c      = length(contractlength);
    temp   = sol.rent;
    disp(['minimal cost of ' num2str(f_k) ' by '])
    for m = 1:Nmonths,
        if sum(temp(:,m)) > 0,
            rent = find(temp(:,m));
            disp(['   renting ' num2str(temp(rent,m)) ' vans for ' ...
                num2str(contractlength(rent)) ' months in month ' num2str(m)])
        end
    end
end

% MODIFICATION LOG
%
% 051021 med   Created.
% 060112 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
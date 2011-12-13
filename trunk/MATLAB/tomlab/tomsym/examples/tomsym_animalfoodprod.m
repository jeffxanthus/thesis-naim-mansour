%% Animal Food Production
%
%% Problem description
% The company CowFood produces food for farm animals that is sold in
% two forms: powder and granules. The raw materials used for the
% production of the food are: oat, maize and molasses. The raw
% materials (with the exception of molasses) first need to be ground,
% and then all raw materials that will form a product are blended. In
% the last step of the production process the product mix is either
% transformed to granules or sieved to obtain food in the form of
% powder.
%
%                    Molasses -+
%                              |
%                              V       +---> Granulating --> Granules
%  Oat ----+                           |
%          +-> Grinding --> Blending --+
%  Maize --+                           |
%                                      +---> Sieving ------> Powder
%  Animal food production process
%
% Every food product needs to fulfill certain nutritional
% requirements. The percentages of proteins, lipids and fibers
% contained in the raw materials and the required percentages in the
% final products are listed in the table below.
%
% Contents of nutritional components in percent
%
%  +-----------------+--------+-------+-----+
%  |Raw material     |Proteins|Lipids |Fiber|
%  +-----------------+--------+-------+-----+
%  |Oat              | 13.6   |  7.1  | 7.0 |
%  |Maize            |  4.1   |  2.4  | 3.7 |
%  |Molasses         |  5.0   |  0.3  | 25  |
%  |Required contents| >=9.5  |  >=2  | <=6 |
%  +-----------------+--------+-------+-----+
%
% There are limits on the availability of raw materials. The table
% below displays the amount of raw material that is available every
% day and the respective prices.
%
% Raw material availabilities and prices
%
%  +------------+----------------------+------------+
%  |Raw material|Available amount in kg|Cost in $/kg|
%  +------------+----------------------+------------+
%  |Oat         | 11900                |   0.13     |
%  |Maize       | 23500                |   0.17     |
%  |Molasses    |   750                |   0.12     |
%  +------------+----------------------+------------+
%
% The cost of the different production steps are given in the
% following table.
%
% Production costs in $/kg
%
%  +--------+--------+-----------+-------+
%  |Grinding|Blending|Granulating|Sieving|
%  +--------+--------+-----------+-------+
%  |  0.25  | 0.05   |    0.42   |  0.17 |
%  +--------+--------+-----------+-------+
%
% With a daily demand of nine tonnes of granules and twelve tonnes of
% powder, which quantities of raw materials are required and how
% should they be blended to minimize the total cost?
%
%% Variables
%
%  compsize                   Amount of different food to produce
%  mincomp                    Minimal content of nutrients
%  maxcomp                    Maximal content of nutrients
%  rawcompmat                 Content of nutrients in raw sources
%  rawavail                   Available raw sources
%  rawcost                    Cost of raw material
%  prodcostsraw               Cost to process raw material
%  prodcostsprod              Cost to produce products
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
compsize    = [9000;12000];
rawcompmat  = [[13.6;4.1;5],...
    [7.1;2.4;.3],...
    [7;3.7;25]];
rawavail     = [11900;23500;750];
rawcost      = [.13;.17;.12];

grindcost = 0.25;
blendcost = 0.05;
grancost = 0.42;
sievcost = 0.17;
use = tom('use',length(rawcost),length(compsize));
objective = sum(sum(repmat(rawcost,1,2).*use)) + ...
    sum(sum(grindcost*use(1:end-1,:))) + ...
    sum(sum(blendcost*use)) + ...
    sum(sum(grancost*use(:,1))) + ...
    sum(sum(sievcost*use(:,2)));

prodcons = {sum(use,1) == compsize'};
reqcons1 = {rawcompmat(:,1)'*use >= 9.5*compsize'};
reqcons2 = {rawcompmat(:,2)'*use >= 2*compsize'};
reqcons3 = {rawcompmat(:,3)'*use <= 6*compsize'};
availcons = {sum(use,2) <= rawavail};
bnds = {0 <= use};
constraints = {prodcons, reqcons1, reqcons2, reqcons3, availcons, bnds};
options = struct;
options.solver = 'cplex';
options.name   = 'Animal Food Production';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    raws  = length(rawavail);
    prods = length(compsize);

    for prd = 1:prods,
        disp(['produce ' num2str(compsize(prd)) ...
            ' of product ' num2str(prd) ','])
        disp('   and use the following ingredients:')

        for raw = 1:raws,
            disp(['   ' num2str(sol.use(raw,prd)) ...
                ' units of ingredient ' num2str(raw) ])
        end
    end
end

% MODIFICATION LOG
%
% 051007 med   Created.
% 060110 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
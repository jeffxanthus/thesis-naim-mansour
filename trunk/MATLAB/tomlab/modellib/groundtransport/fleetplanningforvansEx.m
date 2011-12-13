% function Result = fleetplanningforvansEx(PriLev)
%
% Creates a TOMLAB MIP problem for fleet planning for vans
%
% FLEET PLANNING FOR VANS
%
% A chain of department stores uses a fleet of vans rented from
% different rental agencies. For the next six months period it has
% forecast the following needs for vans (table below):
%
% Requirements for vans for six months
%
% +---+---+---+---+---+---+
% |Jan|Feb|Mar|Apr|May|Jun|
% +---+---+---+---+---+---+
% |430|410|440|390|425|450|
% +---+---+---+---+---+---+
%
% At the 1st January, the chain has 200 vans, for which the rental
% period terminates at the end of February.
%
% To satisfy its needs, the chain has a choice among three types of
% contracts that may start the first day of every month: 3-months
% contracts for a total cost of $1700 per van, 4-months contracts at
% $2200 per van, and 5-months contracts at $2600 per van. How many
% contracts of the different types need to be started every month in
% order to satisfy the company’s needs at the least cost and to have
% no remaining vans rented after the end of June?
%
% Running contracts in month 5 (May)    .........
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
% VARIABLES
%
% demand                     Demand of vans per month
% initialsupply              Vans per month
% contractlength             Contracts available
% contractcost               Cost of contracts
%                             
% RESULTS
%
% To interpret the result, use PriLev > 1, for example:
% Result = fleetplanningforvansEx(2);
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
% Written Oct 21, 2005.   Last modified Oct 21, 2005.

function Result = fleetplanningforvansEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand         = [ 430; 410; 440; 390; 425; 450];

initialsupply  = [ 200; 200;   0;   0;   0;   0];

contractlength = [ 3; 4; 5];

contractcost   = [ 1700; 2200; 2600];

Prob = fleetplanningforvans(demand, initialsupply, contractlength, contractcost);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   months  = length(demand);
   c       = length(contractlength);
   temp   = reshape(Result.x_k,c,months);
   disp(['minimal cost of ' num2str(Result.f_k) ' by '])
   for m = 1:months,
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

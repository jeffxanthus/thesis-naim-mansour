% function Result = familybudgetEx(PriLev)
%
% Creates a TOMLAB MIP problem for family budget
%
% FAMILY BUDGET
%
% The mother of a family wishes to use her son’s computer. On the
% internet she has found optimization software and would now like to
% formulate a mathematical model to help plan their annual budget.
% She has prepared a list of the monthly expenses and receipts of
% money. Every month the expenses for living are $ 550. The monthly
% rent for the apartment is $ 630. She also budgets $ 135 for
% telephone every two months, $ 850 for gas and electricity bills
% every six months, $ 340 per month for the car and $ 100 of tax
% every four months. Her receipts of money are a monthly payment of
% $ 150 of state allowance for families with dependent children and a
% net salary of $ 1900 per month. She knows that she pays at least
% $ 165 for leisure every month (subscription to the swimming pool
% for the older children, football club for the youngest, gym for
% herself) but she would like to spend more (restaurant, cinema,
% holidays). How should the budget be balanced all through the year
% to maximize the money available for leisure?
%
% VARIABLES
%
% costs                      Various costs
% costfreq                   Frequency of the costs
% income                     Various incomes
% incomefreq                 Frequency of the income
% minhobby                   Mimial cost for hobbies
%
% RESULTS
%
% For an interpretation of the results, run:
% Result  = familybudgetEx(2);
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
% Result       Result structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Result = familybudgetEx(PriLev)

if nargin < 1
   PriLev = 1;
end

costs         = [550 630 135 850 340 100]';
costfreq      = [  1   1   2   6   1   4]';
income        = [150 1900]';
incomefreq    = [  1    1]';
minhobby      = 165;

Prob = familybudget(costs, costfreq, income, incomefreq, minhobby);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   m       = 12 ; % number of months
   leisure = Result.x_k(1:m);
   save    = Result.x_k(m+1:2*m);
   disp('Money for leisure and for savings per month')
   for i = 1:m,
      disp(['month ' num2str(i) ', leisure ' num2str(leisure(i))  ', save ' num2str(save(i))])
   end
   disp(['total leisure this year: ' num2str(sum(leisure))])
   disp(['total to save this year: ' num2str(sum(save))   ])
end

% MODIFICATION LOG
%
% 051201 med   Created.
% 060116 per   Added documentation.
% 060117 per   Minor change of documentation.
% 060125 per   Moved disp to end

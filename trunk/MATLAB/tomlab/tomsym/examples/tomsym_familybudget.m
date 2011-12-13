%% Family Budget
%
%% Problem description
% The mother of a family wishes to use her son's computer. On the
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
%% Variables
%
%  costs                      Various costs
%  costfreq                   Frequency of the costs
%  income                     Various incomes
%  incomefreq                 Frequency of the income
%  minhobby                   Minimal cost for hobbies
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
costs      = [550 630 135 850 340 100]';
costfreq   = [  1   1   2   6   1   4]';
income     = [150 1900]';
incomefreq = [  1    1]';
minhobby   = 165;

m = 12; % Hobby, months

hobby = tom('hobby',m,1);
save  = tom('save',m,1);

% No variables are integer.
bnds = {hobby >= minhobby, save >= 0};

% Months constraints
incomefreqcount = incomefreq;
costfreqcount = costfreq;
uppr = zeros(m,1);
value = 0;
for i=1:m
    for j=1:length(incomefreq)
        if incomefreq(j) == i
            value = value + income(j);
            incomefreq(j) = incomefreq(j) + incomefreqcount(j);
        end
    end
    for j=1:length(costfreq)
        if costfreq(j) == i
            value = value - costs(j);
            costfreq(j) = costfreq(j) + costfreqcount(j);
        end
    end
    uppr(i,1) = value;
    value = 0;
end

con1 = {hobby(1) + save(1) <= uppr(1)};
con2 = {hobby(2:end) + save(2:end) - save(1:end-1) <= uppr(2:end)};

% Objective
objective = -sum(hobby);

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Family Budget';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    m       = 12 ; % number of months
    leisure = sol.hobby;
    save    = sol.save;
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
% 090308 med   Converted to tomSym
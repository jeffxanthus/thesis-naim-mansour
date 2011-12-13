%% Choice of Loans
%
%% Problem description
% Mr Chic, director of a chain of shops selling clothes, wishes to
% open three new shops: one in London, one in Munich, and one in
% Rome. To open a new shop costs respectively $ 2.5 million, $ 1
% million and $ 1.7 million. To finance his project, he calls at
% three different banks.
%
% Rates offered by the banks for the different projects
%
%  +------+------+------+----+
%  |      |London|Munich|Rome|
%  +------+------+------+----+
%  |Bank 1| 5.0% | 6.5% |6.1%|
%  |Bank 2| 5.2% | 6.2% |6.2%|
%  |Bank 3| 5.5% | 5.8% |6.5%|
%  +------+------+------+----+
%
% Depending on the location of the shops and the evaluation of the
% related risks, each bank decides to finance at most $ 3 million
% over 8 years with different interest rates for the shops. Determine
% the amount to borrow from each bank for financing each shop in
% order to minimize the total expenses of Mr Chic.
%
%% Variables
%
%  rates                      The rates
%  costs                      Cost per site
%  maxloan                    Max loan per bank
%  loanlength                 Length of the loan in years
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
rates = [ 5.0 6.5 6.1 ;...
    5.2 6.2 6.2 ;...
    5.5 5.8 6.5 ];

costs       = [2.5 1.0 1.7]'*1e6;
maxloan     = 3e6;
loanlength  = 8;

b = size(rates,1);
s = size(rates,2);

borrow = tom('borrow',b,s);

% No variables are binary.
bnds = {0 <= borrow};

% Shop constraints
con1 = {sum(borrow,1)' == costs};

% Bank constraints
con2 = {sum(borrow,2) <= maxloan};

% Objective
rates = rates/100;
objective = sum(sum(borrow.*(rates./(1-(1+rates).^(-loanlength)))));

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Depot Location';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    temp   = sol.borrow;
    for i = 1:s,
        site = temp(:,i);
        idx  = find(site);
        disp(['To finance shop at site ' num2str(i)])
        for j = 1:length(idx),
            disp(['  take a loan of ' num2str(site(idx(j))) ...
                ' in bank ' num2str(idx(j))])
        end
    end
end

% MODIFICATION LOG
%
% 051130 med   Created.
% 060116 per   Added documentation.
% 090308 med   Converted to tomSym
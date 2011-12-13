%% Portfolio Selection
%
%% Problem description
% A consultant in finance has to choose for one of his wealthy female
% clients a certain number of shares in which to invest. She wishes
% to invest $ 100,000 in 6 different shares. The consultant estimates
% for her the return on investment that she may expect for a period
% of six months. The following table gives for each share its country
% of origin, the category (T: technology, N: non-technology) and the
% expected return on investment (ROI). The client specifies certain
% constraints. She wishes to invest at least $ 5,000 and at most
% $ 40,000 into any share. She further wishes to invest half of her
% capital in European shares and at most 30% in technology. How
% should the capital be divided among the shares to obtain the
% highest expected return on investment?
%
% List of shares
%
%  +--+-------+--------+------------+
%  |Nr|Origin |Category|Expected ROI|
%  +--+-------+--------+------------+
%  | 1|Japan  |   T    |   5.3%     |
%  | 2|UK     |   T    |   6.2%     |
%  | 3|France |   T    |   5.1%     |
%  | 4|USA    |   N    |   4.9%     |
%  | 5|Germany|   N    |   6.5%     |
%  | 6|France |   N    |   3.4%     |
%  +--+-------+--------+------------+
%
%% Variables
%
%  budget                     Budget
%  mininvest                  Minimal investment
%  maxinvest                  Maximal investment
%  catinvest1min              Minimal investment in category One (N)
%  idx1cat                    Index of category One
%  catinvest2max              Maximal investment in category Two (T)
%  idx2cat                    Index of category Two
%  returns                    Expected ROI
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
budget        = 100000;
mininvest     = 5000;
maxinvest     = 40000;
catinvest1min = 0.5;
idx1cat       = [0 1 1 0 1 1]';

catinvest2max = 0.3;
idx2cat       = [1 1 1 0 0 0]';

returns       = [5.3 6.2 5.1 4.9 6.5 3.4]';

n   = length(returns);
buy = tom('buy',n,1);

% No variables are binary.
bnds = {mininvest <= buy <= maxinvest};

% Budget constraints
con1 = {sum(buy.*idx2cat) <= catinvest2max*budget};
con2 = {sum(buy.*idx1cat) >= catinvest1min*budget};
con3 = {sum(buy) == budget};

% Objective
objective = -sum(returns.*buy/100);

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Portfolio Selection';
Prob = sym2prob('lp',objective,constraints,[],options);
Prob.MIP.SC = 1:Prob.N; % All variables are semi-continuous

PriLev = 1;
Result = tomRun('cplex', Prob, PriLev);
sol = getSolution(Result);

if PriLev > 0
    invest = sol.buy;
    for i = 1:length(invest),
        if invest(i) ~= 0,
            disp(['invest $ ' num2str(invest(i)) ' in share ' num2str(i)])
        end
    end
end

% MODIFICATION LOG
%
% 051201 med   Created.
% 060116 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
%% Financing an Early Retirement Scheme
%
%% Problem description
% The National Agricultural Bank (NAB) decides to establish an early
% retirement scheme for fifteen employees who are taking early
% retirement. These employees will retire over a period of seven
% years starting from the next year. To finance this early retirement
% scheme, the bank decides to invest in bonds for this seven-year
% period. The necessary amounts to cover the pre-retirement leavers
% are given in the following table; they have to be paid at the
% beginning of every year.
%
% Amounts required every year
%
%  +------------------+----+---+---+---+---+----+---+
%  |Year              |   1|  2|  3|  4|  5|   6|  7|
%  +------------------+----+---+---+---+---+----+---+
%  |Amount (in 1000 $)|1000|600|640|480|760|1020|950|
%  +------------------+----+---+---+---+---+----+---+
%
% For the investments, the bank decides to take three different types
% of bonds, SNCF bonds, Fujitsu bonds, and Treasury bonds. The money
% that is not invested in these bonds is placed as savings with a
% guaranteed yield of 3.2%. The table below lists all information
% concerning the yield and the durations of the bonds and also the
% value of a bond. In this type of bond, it is only possible to buy
% an integer number of bonds, and the invested capital remains locked
% in for the total duration of the bond. Every year, only the
% interest on the capital is distributed. The person in charge of the
% retirement plan decides to buy bonds at the beginning of the first
% year, but not in the following years. How should he organize the
% investments in order to spend the least amount of money to cover
% the projected retirement plan?
%
% Information about loans
%
%  +--------+----------------------+---------+--------+
%  |Loan    |Value of bonds (in k$)| Interest|Duration|
%  +--------+----------------------+---------+--------+
%  |SNCF    |     1.0              |   7.0%  | 5 years|
%  |Fujitsu |     0.8              |   7.0%  | 4 years|
%  |Treasury|     0.5              |   6.5%  | 6 years|
%  +--------+----------------------+---------+--------+
%
%% Variables
%
%  amounts                    Retirement costs
%  baseinterest               Interest from money not invested in bonds
%  values                     Values of a bond
%  interest                   Interest from bonds
%  duration                   Years to place the money
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
amounts       = [1000 600 640 480 760 1020 950]'*1000;
baseinterest  = 3.2;
values        = [1000 800 500]';
interest      = [7.0 7.0 6.5]';
duration      = [5.0 4.0 6.0]';

b  = length(values);    % Bonds to buy
t  = length(amounts)-1; % Invest occasions

buy = tom('buy',b,1,'int');
invest = tom('invest',t,1);
capital = tom('capital',1,1);

% Some variables are integer.
bnds = {capital >= 0, invest >= 0, buy >= 0};

% First year constr
con1 = {capital - sum(values.*buy) - invest(1) == amounts(1)};

% 2nd - 4th year constr.
con2 = cell(3,1);
for i=1:3
    con2{i} = {values'.*interest'/100*buy + ...
        invest(i)*(1+baseinterest/100) - invest(i+1) == ...
        amounts(i+1)};
end

% 5th - 6th year constr.
con3 = cell(2,1);
for i=1:2
    idx1 = find(duration == i+3);
    idx2 = find(duration >= i+4);
    con3{i} = {values(idx1)'.*(1+interest(idx1)'/100)*buy(idx1) + ...
        values(idx2)'.*interest(idx2)'/100*buy(idx2) + ...
        (1+baseinterest/100)*invest(i+3) - invest(i+4) == ...
        amounts(i+4)};
end

% 7th year constr.
idx = find(duration == 6);
con4 = {values(idx)'.*(1+interest(idx)'/100)*buy(idx) + ...
    (1+baseinterest/100)*invest(6) == amounts(7)};

% Objective
objective = capital;

constraints = {bnds, con1, con2, con3, con4};
options = struct;
options.solver = 'cplex';
options.name   = 'Financing an Early Retirement Scheme';
Prob = sym2prob('mip',objective,constraints,[],options);
Prob.MIP.cpxControl.EPGAP = 1e-10;

PriLev = 1;
Result = tomRun('cplex', Prob, PriLev);
sol = getSolution(Result);

if PriLev > 0
    buy_bonds = sol.buy;
    buy_other = sol.invest;

    for i = 1:length(buy_bonds),
        if buy_bonds(i) ~= 0,
            disp(['buy ' num2str(buy_bonds(i)) ' bonds of type ' ...
                num2str(i)])
        end
    end

    for i = 1:length(buy_other),
        if buy_other(i) ~= 0,
            disp(['buy "other things" for $ ' ...
                num2str(buy_other(i)) ' year ' num2str(i)])
        end
    end
end

% MODIFICATION LOG
%
% 051201 med   Created.
% 060116 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
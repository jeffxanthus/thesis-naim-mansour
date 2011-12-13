%% Mean Variance Portfolio Selection
%
%% Problem description
% An investor wishes to invest a certain amount of money. He is
% evaluating four different securities (assets) for his investment.
% The securities are US Treasury Bills (‘T-bills’), a computer
% hardware company, a computer software company, and a high-risk
% investment in a theater production. He estimates the mean yield on
% each dollar he invests in each of the securities, and also adopts
% the Markowitz idea of getting estimates of the variance/covariance
% matrix of estimated returns on the securities. (For example,
% hardware and software company worths tend to move together, but are
% oppositely correlated with the success of theatrical production, as
% people go to the theater more when they have become bored with
% playing with their new computers and computer games.) The return on
% theatrical productions are highly variable, whereas the T-bill
% yield is certain. The estimated returns and the variance/covariance
% matrix are given in the table below.
%
% Estimated returns and variance/covariance matrix
%
%  +----------------+--------+--------+--------+-------+
%  |                |Hardware|Software|Show-biz|T-bills|
%  +----------------+--------+--------+--------+-------+
%  |Estimated return|   8    |   9    |    12  |   7   |
%  +----------------+--------+--------+--------+-------+
%  |Hardware        |   4    |   3    |    -1  |   0   |
%  |Software        |   3    |   6    |     1  |   0   |
%  |Show-biz        |  -1    |   1    |    10  |   0   |
%  |T-bills         |   0    |   0    |     0  |   0   |
%  +----------------+--------+--------+--------+-------+
%
% Question 1: Which investment strategy should the investor adopt to
% minimize the variance subject to getting some specified minimum
% target yield?
%
% Question 2: Which is the least variance investment strategy if the
% investor wants to choose at most two different securities (again
% subject to getting some specified minimum target yield)?
%
%% Variables
%
%  estreturn                  Estimated return of securities
%  covmat                     The variance/covariance matrix
%  target                     Minimum target yield
%  maxassets                  Number of securities
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
estreturn = [8 9 12 7]';
covmat    = [ 4  3 -1  0;...
    3  6  1  0;...
    -1  1 10  0;...
    0  0  0  0];
target    = 8.5;

n = length(estreturn);
buy = tom('buy',n,1);

% Bounds
bnds = {0 <= buy <= 1};

% Cost constraints
con1 = {estreturn'*buy >= target};
con2 = {sum(buy) == 1};

% Objective
objective = 0.5*buy'*covmat*buy;

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Mean Variance Portfolio Selection';
sol1 = ezsolve(objective,constraints,[],options);

buybin = tom('buybin',n,1,'int');
bnds2 = {0 <= buybin <= 1};

maxassets = 2;
con3 = {sum(buybin) <= maxassets};
con4 = {buy - buybin <= 0};

constraints = {bnds, con1, con2, con3, con4};
sol2 = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    names  = ['Hardware';
        'Software';
        'Show-biz';
        'T-bills '];
    disp('Answer to Question 1:')
    for i = 1:length(sol1.buy),
        disp(['   invest ' num2str(sol1.buy(i)*100) ...
            '% of the capital in ' names(i,:) ])
    end
    disp('Answer to Question 2 (limited number of securities):')
    x = sol2.buy;
    x(find(x < 1e-6)) = 0;
    for i = 1:length(x),
        if x(i) ~= 0,
            disp(['   invest ' num2str(x(i)*100) ...
                '% of the capital in ' names(i,:)  ])
        end
    end
end

% MODIFICATION LOG
%
% 051202 med   Created.
% 060117 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
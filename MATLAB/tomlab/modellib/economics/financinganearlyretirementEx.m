% function Result = financinganearlyretirementEx(PriLev)
%
% Creates a TOMLAB MIP problem for financing an early retirement scheme
%
% FINANCING AN EARLY RETIREMENT SCHEME
%
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
% +------------------+----+---+---+---+---+----+---+
% |Year              |   1|  2|  3|  4|  5|   6|  7|
% +------------------+----+---+---+---+---+----+---+
% |Amount (in 1000 $)|1000|600|640|480|760|1020|950|
% +------------------+----+---+---+---+---+----+---+
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
% +--------+----------------------+---------+--------+
% |Loan    |Value of bonds (in k$)| Interest|Duration|
% +--------+----------------------+---------+--------+
% |SNCF    |     1.0              |   7.0%  | 5 years|
% |Fujitsu |     0.8              |   7.0%  | 4 years|
% |Treasury|     0.5              |   6.5%  | 6 years|
% +--------+----------------------+---------+--------+
%
% VARIABLES
%
% amounts                    Retirement costs
% baseinterest               Interest from money not invested in bonds
% values                     Values of a bond
% interest                   Interest from bonds
% duration                   Years to place the money
%
% RESULTS
%
% For an interpretation of the results, try the following:
% Result    = financinganearlyretirementEx(2);
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
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Result = financinganearlyretirementEx(PriLev)

if nargin < 1
   PriLev = 1;
end

amounts       = [1000 600 640 480 760 1020 950]'*1000;
baseinterest  = 3.2;
values        = [1000 800 500]';
interest      = [7.0 7.0 6.5]';
duration      = [5.0 4.0 6.0]';

Prob = financinganearlyretirement(amounts, baseinterest, values, interest, duration);
Prob.MIP.cpxControl.EPGAP = 1e-10;
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   b         = 3; % number of bonds
   y         = 7; % number of years
   buy_bonds = Result.x_k(1:b);
   buy_other = Result.x_k(b+1:b+1+y-1);
   
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

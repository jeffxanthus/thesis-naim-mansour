% function Result = choiceofloansEx(PriLev)
%
% Creates a TOMLAB MIP problem for choice of loans
%
% CHOICE OF LOANS
%
% Mr Chic, director of a chain of shops selling clothes, wishes to
% open three new shops: one in London, one in Munich, and one in
% Rome. To open a new shop costs respectively $ 2.5 million, $ 1
% million and $ 1.7 million. To finance his project, he calls at
% three different banks. 
%
% Rates offered by the banks for the different projects
%
% +------+------+------+----+ 
% |      |London|Munich|Rome|
% +------+------+------+----+ 
% |Bank 1| 5.0% | 6.5% |6.1%|
% |Bank 2| 5.2% | 6.2% |6.2%|
% |Bank 3| 5.5% | 5.8% |6.5%|
% +------+------+------+----+ 
%
% Depending on the location of the shops and the evaluation of the
% related risks, each bank decides to finance at most $ 3 million
% over 8 years with different interest rates for the shops. Determine
% the amount to borrow from each bank for financing each shop in
% order to minimize the total expenses of Mr Chic.
%
% VARIABLES
%
% rates                      The rates
% costs                      Cost per site
% maxloan                    Max loan per bank
% loanlength                 Length of the loan in years
%
% RESULTS
%
% for an interpretation of the results, let PriLev > 1:
% Result = choiceofloansEx(2);
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
% Written Nov 30, 2005.   Last modified Nov 30, 2005.

function Result = choiceofloansEx(PriLev)

if nargin < 1
   PriLev = 1;
end

rates       = [ 5.0 6.5 6.1 ;...
      5.2 6.2 6.2 ;...
      5.5 5.8 6.5 ];

costs       = [2.5 1.0 1.7]'*1e6;
maxloan     = 3e6;
loanlength  = 8;

Prob = choiceofloans(rates, costs, maxloan, loanlength);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   [b,s] = size(rates); % banks/sites
   temp   = Result.x_k;
   temp   = reshape(temp,b,s)';
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

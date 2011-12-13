% function Result = meanvarianceportfolioselecEx(PriLev)
%
% Creates a TOMLAB MIQP problem for mean variance portfolio selection
%
% MEAN VARIANCE PORTFOLIO SELECTION
% 
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
% +----------------+--------+--------+--------+-------+
% |                |Hardware|Software|Show-biz|T-bills|
% +----------------+--------+--------+--------+-------+
% |Estimated return|   8    |   9    |    12  |   7   |
% +----------------+--------+--------+--------+-------+
% |Hardware        |   4    |   3    |    -1  |   0   |
% |Software        |   3    |   6    |     1  |   0   |
% |Show-biz        |  -1    |   1    |    10  |   0   |
% |T-bills         |   0    |   0    |     0  |   0   |
% +----------------+--------+--------+--------+-------+
%
% Question 1: Which investment strategy should the investor adopt to
% minimize the variance subject to getting some specified minimum
% target yield?
%
% Question 2: Which is the least variance investment strategy if the
% investor wants to choose at most two different securities (again
% subject to getting some specified minimum target yield)?
%
% VARIABLES
%
% estreturn                  Estimated return of securities
% covmat                     The variance/covariance matrix
% target                     Minimum target yield
% maxassets                  Number of securities
%
% RESULTS
%
% For an interpretation of the results, try the following:
% [Result1, Result2] = meanvarianceportfolioselecEx(2);
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
% Written Dec 2, 2005.   Last modified Dec 2, 2005.

function [Result1, Result2] = meanvarianceportfolioselecEx(PriLev)

if nargin < 1
   PriLev = 1;
end

estreturn     = [8 9 12 7]';
covmat        = [ 4  3 -1  0;...
      3  6  1  0;...
      -1  1 10  0;...
      0  0  0  0];
target        = 8.5;

Prob = meanvarianceportfolioselec(estreturn, covmat, target);
Result1 = tomRun('cplex', Prob, PriLev);

maxassets     = 2;

Prob = meanvarianceportfolioselec(estreturn, covmat, target, maxassets);
Result2 = tomRun('cplex', Prob, PriLev);

if PriLev >1,
   names  = ['Hardware';
      'Software';
      'Show-biz';
      'T-bills '];
   disp('Answer to Question 1:')
   for i = 1:length(Result1.x_k),
      disp(['   invest ' num2str(Result1.x_k(i)*100) ...
            '% of the capital in ' names(i,:) ]) 
   end
   disp('Answer to Question 2 (limited number of securities):')
   x = Result2.x_k(1:length(Result2.x_k)/2);
   x(find(x < 1e-6)) = 0;
   for i = 1:(length(Result2.x_k)/2),
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

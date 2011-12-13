% function Result = planthepersonnelataconstrEx(PriLev)
%
% Creates a TOMLAB MILP problem for planning the personnel at a construction site
%
% PLANNING THE PERSONNEL AT A CONSTRUCTION SITE
%
% Construction workers who erect the metal skeleton of skyscrapers
% are called steel erectors. The following table lists the
% requirements for steel erectors at a construction site during a
% period of six months. Transfers from other sites to this one are
% possible on the first day of every month and cost $100 per person.
% At the end of every month workers may leave to other sites at a
% transfer cost of $160 per person. It is estimated that
% understaffing as well as overstaffing cost $200 per month per post
% (in the case of unoccupied posts the missing hours have to be
% filled through overtime work).
%
% Monthly requirement for steel erectors
%
% +---+---+---+---+---+---+
% |Mar|Apr|May|Jun|Jul|Aug|
% +---+---+---+---+---+---+
% | 4 | 6 | 7 | 4 | 6 | 2 |
% +---+---+---+---+---+---+
%
% Overtime work is limited to 25% of the hours worked normally. Every
% month, at most three workers may arrive at the site. The departure
% to other sites is limited by agreements with labor unions to 1/3 of
% the total personnel of the month. We suppose that three steel
% erectors are already present on site at the end of February, that
% nobody leaves at the end of February and that three workers need to
% remain on-site at the end of August. Which are the number of
% arrivals and departures every month to minimize the total cost?
%
% VARIABLES
%
% transferin                 Cost to transfer in
% transferout                Cost to transfer out
% staffingdevcost            Cost for over or under employment
% overtimemax                Maximum overtime
% maxtransferin              Maximum amount to transfer in
% maxtransferout             Maximum amount to transfer out
% startstaff                 Starting staff
% endstaff                   Staff required at the end of the period
% demands                    Staff required each month
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = planthepersonnelataconstrEx(2);
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
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Result = planthepersonnelataconstrEx(PriLev)

if nargin < 1
   PriLev = 1;
end

transferin      = 100;
transferout     = 160;
staffingdevcost = 200;

overtimemax     = 0.25;
maxtransferin   = 3;
maxtransferout  = 1/3;

startstaff      = 3;
endstaff        = 3;

demands         = [4 6 7 4 6 2]';

Prob = planthepersonnelataconstr(transferin, transferout,...
   staffingdevcost, overtimemax, maxtransferin, maxtransferout, ...
   startstaff, endstaff, demands);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   months = ['Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'];
   temp   = reshape(Result.x_k,size(months,1),5);
   for m = 1:size(months,1),
      disp(['-- ' months(m,:) ' --' ])
      disp(['  ' num2str(temp(m,1)) ' worker(s) on site']) 
      disp(['  ' num2str(temp(m,2)) ' worker(s) have arrived']) 
      disp(['  ' num2str(temp(m,3)) ' worker(s) will leave']) 
      disp(['  ' num2str(temp(m,4)) ' worker(s) too many']) 
      disp(['  ' num2str(temp(m,5)) ' worker(s) too few']) 
   end
end

% MODIFICATION LOG
%
% 051205 med   Created.
% 060118 per   Added documentation.
% 060126 per   Moved disp to end

% function Result = schedulingnursesEx(PriLev)
%
% Creates a TOMLAB MILP problem for scheduling nurses
%
% SCHEDULING NURSES
%
% Mr. Schedule has been asked to organize the schedule of nurses for
% the Cardiology service at St. Joseph’s hospital. A working day in
% this service is subdivided into twelve periods of two hours. The
% personnel requirement changes from period to period: for instance,
% only a few nurses are required during the night, but the total
% number must be relatively high during the morning to provide
% specific care to the patients. The following table lists the
% personnel requirement for every time period.
%
% Question 1: 
% Determine the minimum number of nurses required to cover all the
% requirements, knowing that a nurse works eight hours per day and
% that she is entitled to a break of two hours after she has worked
% for four hours.
%
% Question 2:
% The service only has 80 nurses, which is not sufficient with the
% given requirements. Mr. Schedule therefore proposes that part of
% the personnel works two additional hours per day. These two
% additional hours follow immediately after the last four hours,
% without any break. Determine the schedule of the nurses in this
% service that minimizes the number of nurses working overtime.
% 
% Personnel requirement per time period
% 
% +------+-------------+------------------------+
% |Number|Time interval|Minimum number of nurses|
% +------+-------------+------------------------+
% |   0  | 00am - 02am |          15            |
% |   1  | 02am - 04am |          15            |
% |   2  | 04am - 06am |          15            |
% |   3  | 06am - 08am |          35            |
% |   4  | 08am - 10am |          40            |
% |   5  | 10am - 12pm |          40            |
% |   6  | 12pm - 02pm |          40            |
% |   7  | 02pm - 04pm |          30            |
% |   8  | 04pm - 06pm |          31            |
% |   9  | 06pm - 08pm |          35            |
% |  10  | 08pm - 10pm |          30            |
% |  11  | 10pm - 12am |          20            |
% +------+-------------+------------------------+
%
%
% VARIABLES
%
% demand                     nurses needed per shift
% hoursperday                hours of work per day
% breakinterval              work this long before break
% breaklength                length of break
% intervallength             length of an interval
% maxdemand                  max number of nurses in Q 2
%
% RESULTS
%
% For an interpretation of the results, run:
% [Result1, Result2] = schedulingnursesEx(2);
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

function [Result1, Result2] = schedulingnursesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

demand         = [15 15 15 35 40 40 40 30 31 35 30 20]';
hoursperday    = 8;
breakinterval  = 4;
breaklength    = 2;
intervallength = 2;
maxdemand      = 80;

flag = 0; % No max on nurses

Prob = schedulingnurses(demand, hoursperday, breakinterval, breaklength, ...
   intervallength, maxdemand, flag);
Result1 = tomRun('cplex', Prob, PriLev);

flag = 1; % Max on nurses
Prob = schedulingnurses(demand, hoursperday, breakinterval, breaklength, ...
   intervallength, maxdemand, flag);
Result2 = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   x1        = Result1.x_k;                      % nurses working
   intervals = length(x1);                       % number of intervals
   work      = [x1 circshift(x1,-1) circshift(x1,-3) circshift(x1,-4)];
   work_1    = sum(work')';                      % total nurses
   n         = Result2.x_k(1:intervals);         % regular shifts
   o         = Result2.x_k(intervals+1:end);     % overtime
   work1     = [n circshift(n,-1) ...            % regular -> 4 intervals 
         circshift(n,-3) circshift(n,-4)];
   work2     = [o circshift(o,-1) ...            % overtime -> 5 intervals
         circshift(o,-3) circshift(o,-4) circshift(o,-5)];
   work_2    = sum(work1')' + sum(work2')';      % total work
   
   disp(['Without overtime we need ' num2str(Result1.f_k) ' nurses.'])
   disp(['   This is the schedule:'])
   for i = 1:intervals,
      disp(['   ' num2str((i-1)*2) '.00: ' num2str(x1(i)) ...
            ' start and ' num2str(work_1(i)) ' are working.' ...
            ' ( demand is ' num2str(demand(i)) ')'  ])
   end
   
   disp(['With overtime we have ' num2str(sum(n)) ...
         ' nurses, and ' num2str(sum(o)) ... 
         ' of them will work overtime.'])
   disp(['   This is the schedule:'])
   for i = 1:intervals,
      disp(['   ' num2str((i-1)*2) '.00: ' num2str(n(i)) ...
            ' start a regular shift, '     num2str(o(i)) ...
            ' start an overtime shift and ' num2str(work_2(i)) ...
            ' are working.' ...
            ' ( demand is ' num2str(demand(i)) ')'  ])
   end
end

% MODIFICATION LOG
%
% 051202 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end

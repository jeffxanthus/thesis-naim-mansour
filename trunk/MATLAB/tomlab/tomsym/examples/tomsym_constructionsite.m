%% Planning the Personnel at a Construction Site
%
%% Problem description
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
%  +---+---+---+---+---+---+
%  |Mar|Apr|May|Jun|Jul|Aug|
%  +---+---+---+---+---+---+
%  | 4 | 6 | 7 | 4 | 6 | 2 |
%  +---+---+---+---+---+---+
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
%% Variables
%
%  transferin                 Cost to transfer in
%  transferout                Cost to transfer out
%  staffingdevcost            Cost for over or under employment
%  overtimemax                Maximum overtime
%  maxtransferin              Maximum amount to transfer in
%  maxtransferout             Maximum amount to transfer out
%  startstaff                 Starting staff
%  endstaff                   Staff required at the end of the period
%  demands                    Staff required each month
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
transferin      = 100;
transferout     = 160;
staffingdevcost = 200;

overtimemax     = 0.25;
maxtransferin   = 3;
maxtransferout  = 1/3;

startstaff      = 3;
endstaff        = 3;

demands         = [4 6 7 4 6 2]';

m = length(demands);
onsite = tom('onsite',m,1,'int');
arrive = tom('arrive',m,1,'int');
leave  = tom('leave',m,1,'int');
over   = tom('over',m,1,'int');
under  = tom('under',m,1,'int');

% Variables are binary
bnds = {onsite >= 0, arrive >= 0, leave >= 0, over >= 0, under >= 0};
bnds = {bnds, arrive <= maxtransferin};

% Start constraint
con1 = {onsite(1) == startstaff + arrive(1)};

% Final constraint
con2 = {endstaff == onsite(end) - leave(end)};

% Intermediate constraints
con3 = {onsite(2:end) == onsite(1:end-1) - ...
    leave(1:end-1) + arrive(2:end)};

% Precense constraints
con4 = {onsite - over + under == demands};

% Overtime constraints
con5 = {under <= onsite*overtimemax};

% Leaving constraints
con6 = {leave <= onsite*maxtransferout};

% Objective
objective = transferin*sum(arrive) + transferout*sum(leave) + ...
    staffingdevcost*sum(over+under);
constraints = {bnds, con1, con2, con3, con4, con5, con6};
options = struct;
options.solver = 'cplex';
options.name   = 'Plan Personnel at a Constr Site';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    monthnames = ['Mar'; 'Apr'; 'May'; 'Jun'; 'Jul'; 'Aug'];
    for m = 1:size(monthnames,1),
        disp(['-- ' monthnames(m,:) ' --' ])
        disp(['  ' num2str(sol.onsite(m)) ' worker(s) on site'])
        disp(['  ' num2str(sol.arrive(m)) ' worker(s) have arrived'])
        disp(['  ' num2str(sol.leave(m)) ' worker(s) will leave'])
        disp(['  ' num2str(sol.over(m)) ' worker(s) too many'])
        disp(['  ' num2str(sol.under(m)) ' worker(s) too few'])
    end
end

% MODIFICATION LOG
%
% 051205 med   Created
% 060118 per   Added documentation
% 060126 per   Moved disp to end
% 090325 med   Converted to tomSym
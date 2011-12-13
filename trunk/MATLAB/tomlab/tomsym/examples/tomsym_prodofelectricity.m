%% Production of Electricity
%
%% Problem description
% Power generators of four different types are available to satisfy the
% daily electricity demands (in megawatts) summarized in the following
% table. We consider a sliding time horizon: the period 10pm-12am of day d
% is followed by the period 0am-6am of day d + 1.
%
% Daily electricity demands (in MW)
%
%  +-------+-------+-------+--------+--------+-------+--------+---------+
%  |Period |0am-6am|6am-9am|9am-12pm|12pm-2pm|2pm-6pm|6pm-10pm|10pm-12am|
%  +-------+-------+-------+--------+--------+-------+--------+---------+
%  |Demand |  12000|  32000|   25000|   36000|  25000|   30000|    18000|
%  +-------+-------+-------+--------+--------+-------+--------+---------+
%
% The power generators of the same type have a maximum capacity and may be
% connected to the network starting from a certain minimal power output.
% They have a start-up cost, a fixed hourly cost for working at minimal
% power, and an hourly cost per additional megawatt for anything beyond
% the minimal output. These data are given in the following table.
%
% Description of power generators
%
%  +---------+-----------+------+--------+--------+------------+--------+
%  |Available|Min. output|Max.  |capacity|Fix cost|Add. MW cost|Start-up|
%  |number   |in MW      |in MW |$/h     |        |$/h         |cost    |
%  +---------+-----------+------+--------+--------+------------+--------+
%  |Type 1   |  10       |  750 |  1750  |  2250  |  2.7       | 5000   |
%  |Type 2   |   4       | 1000 |  1500  |  1800  |  2.2       | 1600   |
%  |Type 3   |   8       | 1200 |  2000  |  3750  |  1.8       | 2400   |
%  |Type 4   |   3       | 1800 |  3500  |  4800  |  3.8       | 1200   |
%  +---------+-----------+------+--------+--------+------------+--------+
%
% A power generator can only be started or stopped at the beginning of a
% time period. As opposed to the start, stopping a power plant does not
% cost anything. At any moment, the working power generators must be able
% to cope with an increase by 20% of the demand forecast. Which power
% generators should be used in every period in order to minimize the total
% daily cost?
%
%% Variables
%
%  demand              MW needed each period
%  available           Generators of each type available
%  mincap              Minimal output of generator
%  maxcap              Maximal output of generator
%  fixcost             Fix cost per hour
%  runningcost         Cost per hour and MW above mincap
%  startcost           Cost for starting a generator
%  periodlengths       Hours in each period
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
demand      = [12;32;25;36;25;30;18]*1000;
available   = [10;4;8;3];
mincap      = [750;1000;1200;1800];
maxcap      = [1750;1500;2000;3500];
fixcost     = [2250;1800;3750;4800];
runningcost = [2.7;2.2;1.8;3.8];
startcost   = [5000;1600;2400;1200];
periodlengths = [6;3;3;2;4;4;2];

n1 = length(demand);
n2 = length(available);

start = tom('start',n1,n2,'int');
wrk  = tom('wrk',n1,n2,'int');
padd  = tom('padd',n1,n2);

% Bounds
bnds1 = {[start;padd] >= 0};
bnds2 = {wrk <= repmat(available',n1,1)};

% Relationship between work and padd
addcap = repmat((maxcap-mincap)',n1,1);
con1 = {padd <= addcap.*wrk};

% Demand in each period constraint
con2 = {sum(repmat(mincap',n1,1).*wrk+padd,2) >= demand};

% 20% additional production possible
con3 = {sum(repmat(maxcap',n1,1).*wrk,2) >= 1.2*demand};

% Relationship between start_pt and work_pt
con4 = {start(1,:) >= wrk(1,:) - wrk(end,:)};
con5 = {start(2:end,:) >= wrk(2:end,:) - wrk(1:end-1,:)};

% Objective
objective = sum(sum(repmat(startcost',n1,1).*start + ...
    repmat(periodlengths,1,n2).*(repmat(fixcost',n1,1).*wrk + ...
    repmat(runningcost',n1,1).*padd)));

constraints = {bnds1, bnds2, con1, con2, con3, con4, con5};
options = struct;
options.solver = 'cplex';
options.name   = 'Production of Electricity';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    intervals = ['00-06';'06-09';'09-12';'12-14';'14-18';'18-22';'22-24'];
    starts    = sol.start;
    running   = sol.wrk;
    extra     = sol.padd;
    [reacs, tim] = size(starts);
    for time = 1:tim,
        disp([' At ' intervals(time,:) '...'])
        for reac = 1:reacs,
            if running(reac,time) > 0,
                disp([   '   there are ' num2str(running(reac,time))  ...
                    ' reactors of type ' num2str(reac) ' running' ])
                if starts(reac,time) > 0,
                    disp(['      (we started ' num2str( starts(reac,time))...
                        ' additional reactors)'])
                end
                if extra(reac,time) > 0,
                    disp(['      (we also produce ' num2str(  extra(reac,time))...
                        ' MW more than the minimal level)'])
                end
            end
        end
    end
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060109 per   Added documentation.
% 060125 per   Moved disp to end
% 090308 med   Converted to tomSym
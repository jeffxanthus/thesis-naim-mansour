%% Scheduling of Telecommunications Via Satellite
%
%% Problem description
% A digital telecommunications system via satellite consists of a
% satellite and a set of stations on earth which serve as interfaces
% with the terrestrial network. With SS-TDMA (Satellite-Switch, Time
% Division Multiple Access) access technology, the satellite divides
% its time among the stations. Consider for example the transmissions
% from four transmitter stations in the US to four receiver stations
% in Europe. The following table gives a possible 4 × 4 traffic
% matrix. TRAFtr is the quantity of data transmitted from station t
% to station r. It is expressed in seconds of transmission duration,
% because all lines have the same constant transmission rate.
%
% Matrix TRAF with its lower bound LB
%
%  +----+--+--+--+--+-------+
%  |TRAF| 1| 2| 3| 4| rowt  |
%  +----+--+--+--+--+-------+
%  |  1 | 0| 7|11|15|   33  |
%  |  2 |15| 8|13| 9|   45  |
%  |  3 |17|12| 6|10|   45  |
%  |  4 | 6|13|15| 4|   38  |
%  +----+--+--+--+--+-------+
%  |colr|38|40|45|38|LB = 45|
%  +----+--+--+--+--+-------+
%
% The satellite has a switch that allows any permutation between the
% four transmitters and the four receivers. The figure below gives an
% example of a permutation connecting the transmitters 1 to 4 to the
% receivers 3, 4, 1, and 2 respectively. These connections allow
% routing a part of the traffic matrix, called a mode. A part of a
% matrix entry transmitted during a mode is a packet of data. A mode
% is thus a 4 × 4 matrix M with at most one non-zero packet per row
% and column and such that Mtr <= TRAFtr for all t, r. To every
% mode corresponds a slice of a schedule as shown in the figure.
%
%  +----+--+--+--+--+         +---------+---------------+
%  |TRAF| 1| 2| 3| 4|         |Stations | Packets       |
%  +----+--+--+--+--+         +---------+---------------+
%  |  1 | 0| 0|11| 0|         | 1 --> 3 |-----11--->    |
%  |  2 | 0| 0| 0| 9|         | 2 --> 4 |------9->      |
%  |  3 |15| 0| 0| 0|         | 3 --> 1 |-----15------->|
%  |  4 | 0|13| 0| 0|         | 4 --> 2 |-----13----->  |
%  +----+--+--+--+--+         +---------+---------------+
%  |colr|38|40|45|38|LB = 45
%  +----+--+--+--+--+
%
% Example of a mode and associated schedule
%
% A valid schedule of transmissions defines a sequence of
% permutations for the on-board switch that routes all the traffic of
% the matrix TRAF. In other words, this boils down to decomposing
% TRAF into a sum of mode matrices. An element of TRAF may be
% fragmented, like TRAF31 that is only partially transmitted in the
% mode represented in the figure above. A fragmented element will
% appear in several packets and several modes of the final
% decomposition. The duration of a mode is the length of its longest
% packet. The objective is to find a schedule with minimal total duration.
%
%% Variables
%
%  traffic                    A matrix describing the traffic
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
PriLev = 1;
if PriLev > 0      % if:   PriLev > 0
    sequence = []; % then: let us catch all x_k
    pmins = [];
end

traffic = [  0   7  11  15;...
    15   8  13   9;...
    17  12   6  10;...
    6  13  15   4];

n1 = size(traffic, 1);
rowsum = sum(traffic,1);
colsum = sum(traffic,2);
LB = max(max([rowsum;colsum']));
q = zeros(n1,n1);
for t=1:n1
    for r=1:n1
        q = min([LB-rowsum(t);LB-colsum(r)]);
        traffic(r,t) = traffic(r,t) + q;
        rowsum(t) = rowsum(t) + q;
        colsum(r) = colsum(r) + q;
    end
end

TQBS = traffic;

flows = tom('flows',n1,n1,'int');
pmin = tom('pmin',1,1);

% Bounds
bnds = {0 <= flows <= 1, pmin >= 0};

while sum(sum(TQBS)) > 0.01
    temp = TQBS > 0;
    temp = spones(temp);
    con1 = {sum(temp.*flows,2) == 1};
    con2 = {sum(temp.*flows,1) == 1};
    con3 = {sum(temp.*TQBS.*flows,2) >= pmin};
    objective = -pmin;
    constraints = {bnds, con1, con2, con3};
    options = struct;
    options.solver = 'cplex';
    options.name   = 'Satellite Scheduling';
    sol = ezsolve(objective,constraints,[],options);

    TQBS = TQBS + sol.flows*(subs(objective,sol));

    if PriLev > 0
        sequence = [sequence sol.flows(:)];
        pmins = [pmins sol.pmin];
    end
end

if PriLev > 0
    for t = 1:size(sequence,2),
        disp(['Transmission ' num2str(t) ' transfers ' num2str(pmins(t)) ' packet(s) of data' ])
        this = reshape(sequence(1:end,t),n1,n1);
        this(find(this < 0.5))  = 0; % remove bad zeros
        for j = 1:n1,
            disp(['   from station ' num2str(j) ' to ' num2str(find(this(j,:))) ])
        end
    end
end

% MODIFICATION LOG
%
% 051122 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end
% 060131 per   cleaned up help text
% 071218 ango  Multiple CPLEX solutions handled gracefully
% 090308 med   Converted to tomSym
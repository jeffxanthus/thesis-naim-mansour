%% Assembly Line Balancing
%
%% Problem description
% An electronics factory produces an amplifier on an assembly line
% with four workstations. An amplifier is assembled in twelve
% operations between which there are certain precedence constraints.
% The following table indicates the duration of every task
% (in minutes) and the list of its immediate predecessors (the
% abbreviation PCB used in this table stands for printed circuit
% board).
%
% The production manager would like to distribute the tasks among the
% four workstations, subject to the precedence constraints, in order
% to balance the line to obtain the shortest possible cycle time,
% that is, the total time required for assembling an amplifier. Every
% task needs to be assigned to a single workstation that has to
% process it without interruption. Every workstation deals with a
% single operation at a time. We talk about cycle time because the
% operations on every workstation will be repeated for every
% amplifier. When an amplifier is finished, the amplifiers at
% stations 1 to 3 advance by one station, and the assembly of a new
% amplifier is started at the first workstation.
%
%  List of tasks and predecessors
%  +----+---------------------------------------+--------+------------+
%  |Task|Description                            |Duration|Predecessors|
%  +----+---------------------------------------+--------+------------+
%  |  1 |Preparing the box                      |    3   |    –       |
%  |  2 |PCB with power module                  |    6   |    1       |
%  |  3 |PCB with pre-amplifier                 |    7   |    1       |
%  |  4 |Filter of the amplifier                |    6   |    2       |
%  |  5 |Push-pull circuit                      |    4   |    2       |
%  |  6 |Connecting the PCBs                    |    8   |   2,3      |
%  |  7 |Integrated circuit of the pre-amplifier|    9   |    3       |
%  |  8 |Adjusting the connections              |   11   |    6       |
%  |  9 |Heat sink of the push-pull             |    2   |  4,5,8     |
%  | 10 |Protective grid                        |   13   |  8,11      |
%  | 11 |Electrostatic protection               |    4   |    7       |
%  | 12 |Putting on the cover                   |    3   |  9,10      |
%  +----+---------------------------------------+--------+------------+
%
%  A precedence graph of the same table with times within ( ).
%
%                   1  (3)
%
%                 /   \
%                /     \
%               /       \
%
%              3 (7)      2 (6)
%
%            /   \      / | \
%           /     \    /  |  \
%          /       \  /   |   \
%
%         7 (9)      6    5    4 (6)
%                     (8)  (4)
%         |          |    |    |
%         |          |    |    |
%         |          |    |    |
%                         |    |
%        11 (4)      8    |    |
%                     (11)|    /
%           \      /  \   |   /
%            \    /    \  |  /
%             \  /      \ | /
%
%              10 (13)    9 (2)
%
%                \       /
%                 \     /
%                  \   /
%
%                    12 (3)
%
%% Variables
%
%  duration                   Time or each task.
%  stations                   Stations available.
%  dependsvec1 and 2          For a task t: If dependsvec1[i] = t and
%                             dependsvec2[i] = d then t depends on d.
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written Oct 7, 2005.   Last modified Dec 14, 2009.

%% Problem setup
duration    = [3;6;7;6;4;8;9;11;2;13;4;3];
stations    = 4;
dependsvec1 = [2;3;4;5;6;6;7;8;9;9;9;10;10;11;12;12];
dependsvec2 = [1;1;2;2;2;3;3;6;4;5;8; 8;11; 7; 9;10];

i = length(duration);
m = stations;
proc  = tom('proc',i,m,'int');
cycle = tom('cycle',1,1);

% All slots are integers
bnds = {0 <= proc <= 1, cycle >= 0};

% Every process has to go to one machine
con1 = {sum(proc,2) == 1};

% Precedence constraint
iter = length(dependsvec1);
con2 = {};
for i=1:iter
    idxi = dependsvec1(i);
    idxj = dependsvec2(i);
    con2{i} = {sum(proc(idxi,:).*(1:4)) >= ...
        sum(proc(idxj,:).*(1:4))};
end

% Cycle constraint
con3 = {duration'*proc <= cycle};

% Objective
objective = cycle;

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Assembly Line Balancing';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    t      = length(duration); % number of tasks
    s      = stations;         % number of workstations
    temp   = sol.proc;   % reshaping distribution
    time   = sol.cycle;  % total time

    for i = 1:s,
        disp(['tasks managed by station ' ...  % filter + disp
            num2str(i) ': ' ...
            num2str(find(temp(:,i))')])
    end

    disp(['total time ' num2str(time)])      % disp the time
end

% MODIFICATION LOG
%
% 051010 med   Created
% 060111 per   Added documentation
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
% function Result = assemblylinebalancingEx(PriLev)
%
% Creates a TOMLAB MIP problem for assembly line balancing
%
% Assembly line balancing
%
% An electronics factory produces an amplifier on an assembly line
% with four workstations. An amplifier is assembled in twelve
% operations between which there are certain precedence constraints.
% The following table indicates the duration of every task
% (in minutes) and the list of its immediate predecessors (the
% abbreviation PCB used in this table stands for ‘printed circuit
% board’). 
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
% List of tasks and predecessors
% +----+---------------------------------------+--------+------------+
% |Task|Description                            |Duration|Predecessors|
% +----+---------------------------------------+--------+------------+
% |  1 |Preparing the box                      |    3   |    –       |
% |  2 |PCB with power module                  |    6   |    1       |
% |  3 |PCB with pre-amplifier                 |    7   |    1       |
% |  4 |Filter of the amplifier                |    6   |    2       |
% |  5 |Push-pull circuit                      |    4   |    2       |
% |  6 |Connecting the PCBs                    |    8   |   2,3      |
% |  7 |Integrated circuit of the pre-amplifier|    9   |    3       |
% |  8 |Adjusting the connections              |   11   |    6       |
% |  9 |Heat sink of the push-pull             |    2   |  4,5,8     |
% | 10 |Protective grid                        |   13   |  8,11      |
% | 11 |Electrostatic protection               |    4   |    7       |
% | 12 |Putting on the cover                   |    3   |  9,10      |
% +----+---------------------------------------+--------+------------+
%
% A precedence graph of the same table with times within ( ).
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
% VARIABLES
%
% duration                   Time or each task.
% stations                   Stations available.
% dependsvec1 and 2          For a task t: If dependsvec1[i] = t and
%                            dependsvec2[i] = d then t depends on d.
%
% RESULTS
%
% To get the results interpreted: 
% Result = assemblylinebalancingEx(2);
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
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 17, 2005.   Last modified Jan 2, 2005.

function Result = assemblylinebalancingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

duration    = [3;6;7;6;4;8;9;11;2;13;4;3];
stations    = 4;
dependsvec1 = [2;3;4;5;6;6;7;8;9;9;9;10;10;11;12;12];
dependsvec2 = [1;1;2;2;2;3;3;6;4;5;8; 8;11; 7; 9;10];

Prob = assemblylinebalancing(stations, duration, dependsvec1, dependsvec2);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   t      = length(duration);               % number of tasks
   s      = stations;                       % number of workstations
   temp   = reshape(Result.x_k(1:t*s),t,s); % reshaping distribution
   time   = Result.x_k(t*s+1);              % total time
   
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
% 060111 per   Added documentation.
% 060126 per   Moved disp to end

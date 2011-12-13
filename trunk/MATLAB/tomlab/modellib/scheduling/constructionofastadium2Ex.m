% function Result = constructionofastadium2Ex(PriLev)
%
% Creates a TOMLAB MILP problem for production of a stadium
%
% CONSTRUCTION OF A STADIUM 2
%  
% A town council wishes to construct a small stadium in order to
% improve the services provided to the people living in the
% district. After the invitation to tender, a local construction
% company is awarded the contract and wishes to complete the task
% within the shortest possible time. All the major tasks are listed
% in the following table. The durations are expressed in weeks. Some
% tasks can only start after the completion of certain other tasks.
% The last two columns of the table refer to question 2 which we
% shall see later.
% 
% Data for stadium construction
%
% +----+--------------------------------+--------+------------+-------+---------------+
% |    |                                |        |            |  Max. | Add. cost per |
% |Task|Description                     |Duration|Predecessors|reduct.|week (in 1000$)|
% +----+--------------------------------+--------+------------+-------+---------------+
% |  1 |Installing the construction site|   2    |    none    |   0   |      –        |
% |  2 |Terracing                       |  16    |    1       |   3   |     30        |
% |  3 |Constructing the foundations    |   9    |    2       |   1   |     26        |
% |  4 |Access roads and other networks |   8    |    2       |   2   |     12        |
% |  5 |Erecting the basement           |  10    |    3       |   2   |     17        |
% |  6 |Main floor                      |   6    |   4,5      |   1   |     15        |
% |  7 |Dividing up the changing rooms  |   2    |    4       |   1   |      8        |
% |  8 |Electrifying the terraces       |   2    |    6       |   0   |      –        |
% |  9 |Constructing the roof           |   9    |   4,6      |   2   |     42        |
% | 10 |Lighting of the stadium         |   5    |    4       |   1   |     21        |
% | 11 |Installing the terraces         |   3    |    6       |   1   |     18        |
% | 12 |Sealing the roof                |   2    |    9       |   0   |      –        |
% | 13 |Finishing the changing rooms    |   1    |    7       |   0   |      –        |
% | 14 |Constructing the ticket office  |   7    |    2       |   2   |     22        |
% | 15 |Secondary access roads          |   4    |  4,14      |   2   |     12        |
% | 16 |Means of signalling             |   3    | 8,11,14    |   1   |      6        |
% | 17 |Lawn and sport accessories      |   9    |    12      |   3   |     16        |
% | 18 |Handing over the building       |   1    |    17      |   0   |      –        |
% +----+--------------------------------+--------+------------+-------+---------------+
%
%
% Precedence graph of construction tasks
%
%
%                                 1
%
%                                 |
%
%                                 2
%
%                              /  |  \
%                             /   |   \
%                            /    |    \
%                           /           \
%                          /      3      \
%                         /               \
%                        /        |        \
%                       /         |         \
%                      /          |          \
%                      
%                   14            5     +----- 4
%                     \                /    /    
%                 /    |          |   /    /  /| \
%                /  +--+----------+---    /  / |  \
%               /  /   |          |   +--+  /  |   \
%                -+    |             /     /   |
%            15        |          6       /    |    7
%                      |                 +     |
%             |        |      /   |   \  |     |    |
%             |        |     /    |    \ |     |    |
%             |        |
%             |        \  11      8      9    10   13
%              \        \
%               \        \   \ /         |     |    |
%                \        \   V          |     |    |
%                 \        \                   |    |
%                  \         16         12     |    |
%                   \                          |    |
%                    \        |          |     |    |
%                     \       |                |   / 
%                      \      |         17     |  /  
%                       \     |                | /   
%                        \    |          |    / /    
%                         \   |              / /  
%                          \  |         18  / /   
%                           \ |            / /    
%                            \ \       /  / /     
%                             \ \     /  / /
%                              \ \   /  / /
%                               \      / /
%                                  F    
%
% Question 1: (see constructionofastadium1Ex.m)
% Which is the earliest possible date of completing the construction?
%
% Question 2: (answered here)
% The town council would like the project to terminate earlier than
% the time announced by the builder (answer to question 1). To obtain
% this, the council is prepared to pay a bonus of $30K for every week
% the work finishes early. The builder needs to employ additional
% workers and rent more equipment to cut down on the total time. In
% the preceding table he has summarized the maximum number of weeks
% he can save per task (column "Max. reduct.") and the associated
% additional cost per week. When will the project be % completed if
% the builder wishes to maximize his profit?
%
% VARIABLES
%
% taskduration               The time it takes to complete each task.
% taskprecedence             Describes the order of the tasks.
%
% RESULTS
%
% For an interpretation of the results, use PriLev > 1, for example:
% Result = constructionofastadium2Ex(2);
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
% Written Oct 10, 2005.   Last modified Jan 2, 2006.

function Result = constructionofastadium2Ex(PriLev)

if nargin < 1
   PriLev = 1;
end

taskduration   = [2;16;9;8;10;6;2;2;9;5;3;2;1;7;4;3;9;1];
taskprecedence = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0;...
      0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0;...
      0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;...
      0  0  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0;...
      0  0  0  0  0  0  0  1  0  0  1  0  0  1  0  0  0  0;...
      0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0;...
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];

Prob = constructionofastadium1(taskduration, taskprecedence);
Prob.PriLevOpt = PriLev-1;
Result = tomRun('cplex', Prob, PriLev-1);

maxreduc       = [0;3;1;2;2;1;1;0;2;1;1;0;0;2;2;1;3;0];
costperweek    = [0;30;26;12;17;15;8;0;42;21;18;0;0;22;12;6;16;0]*1000;
bonus          = 30000;
idx            = sum(taskprecedence,2);

Prob = constructionofastadium2(maxreduc, costperweek, bonus, idx, taskduration, Result);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   tasks = length(maxreduc) + 1; % number of tasks plus one 
   temp  = reshape(Result.x_k,tasks,2)';
   [finished, task] = sort(temp(2,:));
   
   disp('for a best solution')
   for i = 1:tasks,
      if temp(1,task(i)) > 0,
         disp(['   finish task ' num2str(task(i)) ' week ' ...
               num2str(temp(2,task(i))) ' (with a reduction of ' ...
               num2str(temp(1,task(i))) ')'])
      else
         disp(['   finish task ' num2str(task(i)) ' week ' ...
               num2str(temp(2,task(i))) ])
      end
   end
end


% MODIFICATION LOG
%
% 051010 med   Created.
% 060111 per   Added documentation.
% 060126 per   Moved disp to end

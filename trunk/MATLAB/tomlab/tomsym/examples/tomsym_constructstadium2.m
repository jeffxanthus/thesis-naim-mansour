%% Construction of a Stadium 2
%
%% Problem description
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
%  +----+---------------------------+----+-------+-----+--------------+
%  |    |                           |    |       | Max.| Add. cost per|
%  |Task|Description                |Dur.| Pred. | red.| wk (in 1000$)|
%  +----+---------------------------+----+-------+-----+--------------+
%  |  1 |Inst. the constr. site     | 2  | none  |  0  |      –       |
%  |  2 |Terracing                  |16  |  1    |  3  |     30       |
%  |  3 |Constructing foundations   | 9  |  2    |  1  |     26       |
%  |  4 |Access roads networks      | 8  |  2    |  2  |     12       |
%  |  5 |Erecting the basement      |10  |  3    |  2  |     17       |
%  |  6 |Main floor                 | 6  | 4,5   |  1  |     15       |
%  |  7 |Dividing up changing rms   | 2  |  4    |  1  |      8       |
%  |  8 |Electrifying the terraces  | 2  |  6    |  0  |      –       |
%  |  9 |Constructing the roof      | 9  | 4,6   |  2  |     42       |
%  | 10 |Lighting of the stadium    | 5  |  4    |  1  |     21       |
%  | 11 |Installing the terraces    | 3  |  6    |  1  |     18       |
%  | 12 |Sealing the roof           | 2  |  9    |  0  |      –       |
%  | 13 |Finishing the changing rms | 1  |  7    |  0  |      –       |
%  | 14 |Constructing ticket office | 7  |  2    |  2  |     22       |
%  | 15 |Secondary access roads     | 4  | 4,14  |  2  |     12       |
%  | 16 |Means of signalling        | 3  |8,11,14|  1  |      6       |
%  | 17 |Lawn and sport accessories | 9  |  12   |  3  |     16       |
%  | 18 |Handing over the building  | 1  |  17   |  0  |      –       |
%  +----+---------------------------+----+-------+-----+--------------+
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
% Question 1: (see tomsym_constructstadium1.m)
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
%% Variables
%
%  taskduration       The time it takes to complete each task.
%  taskprecedence     Describes the order of the tasks.
%
%% References
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 10, 2005.   Last modified Apr 11, 2009.

tomsym_constructstadium1;

maxreduc    = [0;3;1;2;2;1;1;0;2;1;1;0;0;2;2;1;3;0];
costperweek = [0;30;26;12;17;15;8;0;42;21;18;0;0;22;12;6;16;0]*1000;
bonus       = 30000;

n = length(maxreduc);

save = tom('save',n+1,1,'int');
start = tom('start',n+1,1,'int');

obj1 = Result.f_k;

% Save is first variable, start second
bnds = {0 <= save <= [maxreduc;inf],...
    0 <= start};

count = 1;
con1 = {};
for i=1:n
    idx = find(taskprecedence(i,:) ~= 0);
    for j=1:length(idx)
        con1{count} = {start(idx(j)) + taskduration(i) ...
            - save(i) <= start(i)};
        count = count + 1;
    end
end
con1{count} = {start(end-1) <= start(end)};

idx = find(sum(taskprecedence,2) == 0);
count = 1;
con2 = cell(length(idx),1);
for i=1:length(idx)
    con2{count} = {start(i) >= taskduration(idx(i)) - save(i)};
    count = count + 1;
end

% Start-save constraint
con3 = {start(end) == obj1 - save(end)};

% Objective
objective = -bonus*save(end)+costperweek'*save(1:end-1);

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Construction of a Stadium 2';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    tasks = length(maxreduc) + 1; % number of tasks plus one
    temp  = [sol.save'; sol.start'];
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
% 051010 med   Created
% 060111 per   Added documentation
% 060126 per   Moved disp to end
% 090411 med   Converted to tomSym
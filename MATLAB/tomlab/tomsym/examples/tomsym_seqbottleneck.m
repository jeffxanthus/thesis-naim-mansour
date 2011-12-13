%% Sequencing Jobs on a Bottleneck Machine
%
%% Problem description
% In workshops it frequently happens that a single machine determines
% the throughput of the entire production (for example, a machine of
% which only one is available, the slowest machine in a production
% line, etc.). This machine is called the critical machine or the
% bottleneck. In such a case it is important to schedule the tasks
% that need to be performed by this machine in the best possible way.
% The aim of the problem is to provide a simple model for scheduling
% operations on a single machine and that may be used with different
% objective functions. We shall see here how to minimize the total
% processing time, the average processing time, and the total
% tardiness. A set of tasks (or jobs) is to be processed on a single
% machine. The execution of tasks is non-preemptive (that is, an
% operation may not be interrupted before its completion). For every
% task i its release date and duration are given. For the last
% optimization criterion (total tardiness), a due date (latest
% completion time) is also required to measure the tardiness, that
% is, the amount of time by which the completion of jobs exceeds
% their respective due dates. The following table lists the data for
% our problem.
%
% What is the optimal value for each of the objectives:
%    1 - minimizing the total duration of the schedule (= makespan)?
%    2 - minimizing the mean processing time?
%    3 - minimizing the total tardiness?
%
% Task time windows and durations
%
%  +------------+--+--+--+--+--+--+--+
%  |Job         | 1| 2| 3| 4| 5| 6| 7|
%  +------------+--+--+--+--+--+--+--+
%  |Release date| 2| 5| 4| 0| 0| 8| 9|
%  |Duration    | 5| 6| 8| 4| 2| 4| 2|
%  |Due date    |10|21|15|10| 5|15|22|
%  +------------+--+--+--+--+--+--+--+
%
%% Variables
%
%  releasedate        Release dates of jobs (row 1 in table)
%  duration           Time to perform a job (row 2 in table)
%  duedate            Deadline of each job  (row 3 in table)
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
releasedate = [2 5 4 0 0 8 9]';
duration    = [5 6 8 4 2 4 2]';
duedate     = [10 21 15 10 5 15 22]';

j = 7;
k = 7;

rnk = tom('rnk',j,k,'int');
start = tom('start',k,1,'int');

% All slots are integers
bnds = {0 <= rnk <= 1, start >= 0};

% Assignment constraint one job per rank.
con1 = {sum(rnk,2) == 1};

% Assignment constraint one rank per job.
con2 = {sum(rnk,1) == 1};

% Release constraints.
con3 = {start >= (releasedate'*rnk)'};

% No simultaneous constraints
con4 = {start(2:end) >= start(1:end-1) + (duration'*rnk(:,1:end-1))'};

% Objective 1
objective = start(end) + sum(duration'*rnk(:,end));

constraints = {bnds, con1, con2, con3, con4};
options = struct;
options.solver = 'cplex';
options.name   = 'Sequencing with Bottleneck';
sol1 = ezsolve(objective,constraints,[],options);

% Objective 2
comp = tom('comp',k,1,'int');
bnds2 = {comp >= 0};
con5 = {comp == start + (duration'*rnk)'};
objective = sum(comp);
constraints = {constraints, bnds2, con5};
sol2 = ezsolve(objective,constraints,[],options);

% Objective 3
late = tom('late',k,1,'int');
bnds3 = {late >= 0};
con6 = {late >= comp - (duedate'*rnk)'};
constraints = {constraints, bnds3, con6};
objective = sum(late);
sol3 = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    j         = length(releasedate); % number of jobs
    sequence2 = [sol2.rnk, sol2.start, sol2.comp]; % results from 1 and 2
    sequence3 = [sol3.rnk, sol3.start, sol3.comp, sol3.late]; % results from 3
    sequence2(find(sequence2(1:7*7)<0.1)) = 0; % remove bad zeros
    sequence3(find(sequence3(1:7*7)<0.1)) = 0; % remove bad zeros
    s2 = [];                                   % blank sequence
    s3 = [];                                   % blank sequence
    for t = 1:j,
        s2 = [s2 find(sequence2(t,1:j))];
        s3 = [s3 find(sequence3(t,1:j))];
    end
    disp(['An order to minimize duration (=' ...
        num2str(sequence2(j,j+2)) ') and to minimize mean '...
        'processing time (=' num2str(sequence2(j,j+2)/j) ')'])
    disp(num2str(s2))
    disp(' ')
    tard = sum(sequence3(:,size(sequence3,2)));
    disp(['An order to minimize tardiness (=' num2str(tard) ')' ])
    disp(num2str(s3))
    disp(' ')
end

% MODIFICATION LOG
%
% 060102 med   Created.
% 060111 per   Added documentation.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
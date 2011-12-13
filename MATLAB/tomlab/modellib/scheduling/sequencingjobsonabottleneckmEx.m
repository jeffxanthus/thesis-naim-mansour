% function Result = sequencingjobsonabottleneckmEx(PriLev)
%
% Creates a TOMLAB MIP problem for sequencing jobs on a bottleneck machine
% 
% SEQUENCING JOBS ON A BOTTLENECK MACHINE
%
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
% +------------+--+--+--+--+--+--+--+
% |Job         | 1| 2| 3| 4| 5| 6| 7|
% +------------+--+--+--+--+--+--+--+
% |Release date| 2| 5| 4| 0| 0| 8| 9|
% |Duration    | 5| 6| 8| 4| 2| 4| 2|
% |Due date    |10|21|15|10| 5|15|22|
% +------------+--+--+--+--+--+--+--+
%
% VARIABLES
%
% releasedate                Release dates of jobs (row 1 in table)
% duration                   Time to perform a job (row 2 in table)
% duedate                    Deadline of each job  (row 3 in table)
%
% RESULTS
%
% For an interpretation of the results try the following:
% [Result1, Result2, Result3] = sequencingjobsonabottleneckmEx(2)
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
% Written Jan 2, 2006.   Last modified Jan 2, 2006.

function [Result1,Result2,Result3] = sequencingjobsonabottleneckmEx(PriLev)

if nargin < 1
   PriLev = 1;
end

releasedate = [2 5 4 0 0 8 9]';
duration    = [5 6 8 4 2 4 2]';
duedate     = [10 21 15 10 5 15 22]';

Prob1 = sequencingjobsonabottleneckm(releasedate, duration, duedate, 1);
Result1 = tomRun('cplex', Prob1, PriLev);

Prob2 = sequencingjobsonabottleneckm(releasedate, duration, duedate, 2);
Result2 = tomRun('cplex', Prob2, PriLev);

Prob3 = sequencingjobsonabottleneckm(releasedate, duration, duedate, 3);
Result3 = tomRun('cplex', Prob3, PriLev);

if PriLev > 1,
   j         = length(releasedate); % number of jobs
   sequence2 = reshape(Result2.x_k,j,j+2);    % results from 1 and 2
   sequence3 = reshape(Result3.x_k,j,j+3);    % results from 3
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

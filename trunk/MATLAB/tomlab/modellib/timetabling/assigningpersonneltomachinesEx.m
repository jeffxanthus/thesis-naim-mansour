% function Result = assigningpersonneltomachinesEx(PriLev)
%
% Creates a TOMLAB MILP problem for assigning personnel to machines
%
% ASSIGNING PERSONNEL TO MACHINES
%
% An operator needs to be assigned to each of the six machines in a
% workshop. Six workers have been pre-selected. Everyone has
% undergone a test of her productivity on every machine. The table
% below lists the productivities in pieces per hour. The machines run
% in parallel, that is, the total productivity of the workshop is the
% sum of the productivities of the people assigned to the machines.
%
% Productivity in pieces per hour
%
% +-------+-----------------+ 
% |       |    Machines     |
% +-------+--+--+--+--+--+--+ 
% |Workers| 1| 2| 3| 4| 5| 6|
% +-------+--+--+--+--+--+--+ 
% |   1   |13|24|31|19|40|29|
% |   2   |18|25|30|15|43|22|
% |   3   |20|20|27|25|34|33|
% |   4   |23|26|28|18|37|30|
% |   5   |28|33|34|17|38|20|
% |   6   |19|36|25|27|45|24|
% +-------+--+--+--+--+--+--+ 
%
% The objective is to determine an assignment of workers to machines
% that maximizes the total productivity. We may start by calculating
% a (non-optimal) heuristic solution using the following fairly
% natural method: choose the assignment p -> m with the highest
% productivity, cross out the line p and the column m (since the
% person has been placed and the machine has an operator), and
% restart this process until we have assigned all persons. The
% problem should then be solved to optimality using Mathematical
% Programming. And finally, solve the same problem to optimality,
% but for machines working in series.
%
% VARIABLES
%
% prodmat                    The productivity matrix
%
% RESULTS
%
% For an interpretation of the results, run:
% [Result1, Result2] = assigningpersonneltomachinesEx(2);
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

function [Result1, Result2] = assigningpersonneltomachinesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

prodmat       = [ 13 24 31 19 40 29;...
      18 25 30 15 43 22;...
      20 20 27 25 34 33;...
      23 26 28 18 37 30;...
      28 33 34 17 38 20;...
      19 36 25 27 45 24];

flag = 0; % Parallel

Prob = assigningpersonneltomachines(prodmat, flag);
Result1 = tomRun('cplex', Prob, PriLev);

flag = 1; % Series

Prob = assigningpersonneltomachines(prodmat, flag);
Result2 = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   m   = size(prodmat,1); % number of machines
   w   = m;               % number of workers
   x1  = reshape(Result1.x_k,m,w);
   disp(['Best parallel work (' num2str(-Result1.f_k) ') when '])
   [worker,machine] = find(x1);
   for i = 1:length(worker)
      disp(['   worker '           num2str(worker(i)) ...
            ' operates machine ' num2str(machine(i))])
   end
   x2  = reshape(Result2.x_k(1:m*w),m,w);
   disp(['Best serial work (' num2str(-Result2.f_k) ') when '])
   [worker,machine] = find(x2);
   for i = 1:length(worker)
      disp(['   worker '           num2str(worker(i)) ...
            ' operates machine ' num2str(machine(i))])
   end
end

% MODIFICATION LOG
%
% 051202 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end

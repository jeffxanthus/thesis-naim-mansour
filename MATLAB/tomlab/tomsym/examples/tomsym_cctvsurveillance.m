%% CCTV Surveillance
%
%% Problem description
% In the course of the last few months, the industrial zone of
% Billston has suffered from a series of break-ins and thefts during
% the night. The zone is watched by security men but there are too
% few of them. The town council in charge of security in this zone
% decides to install surveillance cameras to aid the security men
% with their task. These cameras can be directed and pivot through
% 360 degrees. By installing a camera at the intersection of several
% streets, it is possible to survey all adjoining streets. The map in
% the figure below shows the industrial zone with the limits of the
% zone to be covered by closed circuit TV (CCTV) surveillance and the
% 49 possible locations where to install the cameras. What is the
% minimum number of cameras that have to be installed to survey all
% the streets of this zone and where should they be placed?
%
% The industrial zone in Billston
%
%  13 -- 14 -- 18 -- 17    28 -- 29    35 -- 36
% 
%   |     |     |           |           |
%   |     |     |           |           |
%   |
%   |    15 -- 19          26 -- 27    34    48
%   |
%   |  /  |     |           |           |     |
%   | /   |     |           |           |     |
% 
%  12    16 -- 20    24 -- 25 -- 30    33    47 -- 45 -- 46
% 
%   |  /        |        /        |  /  |           |
%   | /         |       /         | /   |           |
% 
%   3 -- 11 -- 21 -- 22          31    37 -- 43 -- 44 -- 49
% 
%   | \                 \         |     |
%   |  \                 \        |     |
%   |
%   |    4 -- 9 -- 10      23 -- 32 -- 38
%   |
%   |    | \                      |     |
%   |    |  \                     |     |
%   |
%   |    6     5                 39 -- 40 -- 41 -- 42
%   |
%   |    | \                      |        /
%   |    |  \                     |       /
%   |                             |      /
%   |    8     7                  |     /
%   |                             |    /
%   |                             |   /
%   |                             |  /
%   |                             | /
% 
%   1 --------------------------- 2
%
%% Variables
%
%  arcs_out/in                These variables describe the network
%                             of streets
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Mar 24, 2009.

%% Problem setup
arcs_in =  [1 1 2 2 3 3 3 3 4 4 4 6 6 9 11 12 12 13 14 14 15 15 16 ...
    17 18 19 20 21 22 22 23 24 25 25 26 26 28 30 31 31 32 32 ...
    33 33 34 35 37 37 38 39 40 41 43 44 44 45 45 47]';

arcs_out = [2 3 39 41 4 11 12 16 5 6 9 7 8 10 21 13 15 14 15 18 16 ...
    19 20 18 19 20 21 22 23 25 32 25 26 30 27 28 29 31 32 33 ...
    38 39 34 37 35 36 38 43 40 40 41 42 44 49 45 46 47 48]';

n = length(arcs_in);  %arcs
arcs = tom('arcs',n,1,'int');

% All variables are binary
bnds = {0 <= arcs <= 1};

% All streets need to be covered
con = {arcs(arcs_in) + arcs(arcs_out) >= 1};

% Objective
objective = sum(arcs);
constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'CCTV Surveillance';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    disp('Put cameras in nodes ')
    disp(num2str(find(sol.arcs)'))
end

% MODIFICATION LOG
%
% 051205 med   Created
% 060118 per   Added documentation
% 060125 per   Moved disp to end
% 090325 med   Converted to tomSym
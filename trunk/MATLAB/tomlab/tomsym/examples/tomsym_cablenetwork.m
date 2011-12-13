%% Construction of a Cabled Network
%
%% Problem description
% A university wishes to connect six terminals located in different
% buildings of its campus. The distances, in meters, between the
% different terminals are given in the table below.
%
% Distances between the different terminals (in meters)
%
%  +----------+---+---+---+---+---+---+
%  |          | T1| T2| T3| T4| T5| T6|
%  +----------+---+---+---+---+---+---+
%  |Terminal 1|  0|120| 92|265|149|194|
%  |Terminal 2|120|  0|141|170| 93|164|
%  |Terminal 3| 92|141|  0|218|103|116|
%  |Terminal 4|265|170|218|  0|110|126|
%  |Terminal 5|149| 93|103|110|  0| 72|
%  |Terminal 6|194|164|116|126| 72|  0|
%  +----------+---+---+---+---+---+---+
%
% These terminals are to be connected via underground cables. We
% suppose the cost of connecting two terminals is proportional to the
% distance between them. Determine the connections to install to
% minimize the total cost.
%
%% Variables
%
%  distances                  a distance matrix
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
distances    = [  0 120  92 265 149 194;...
    120   0 141 170  93 164;...
    92 141   0 218 103 116;...
    265 170 218   0 110 126;...
    149  93 103 110   0  72;...
    194 164 116 126  72   0];

s = size(distances, 1);
t = s;

connected = tom('connected',s,t,'int');
level   = tom('level',s,1);

% Some variables are binary.
bnds = {0 <= connected <= 1, level >= 0, connected(1:s+1:s^2) == 0};

% Terminal connections constraint
con1 = {sum(sum(connected)) <= s-1};

% At least one terminal connected to other
con2 = cell(s*t-s,1);
count = 1;
for i=1:s
    for j=1:t
        if i~=j
            con2{count} = level(j) >= level(i) + 1 - s + s*connected(i,j);
            count = count + 1;
        end
    end
end

% Anti-cycling constraint
con3 = {sum(connected(2:end,:),2) == 1};

% Objective
objective = sum(sum(distances.*connected));

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Construction of a Cabled Network';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    t      = size(distances,1);    % number of terminals
    s      = sum(1:t-1);           % possible connections
    arcs   = [];                   % empty set of arcs
    count  = 1;                    % counter
    % we catch true arcs
    [rows,cols] = find(sol.connected);
    for i = 1:length(rows),
        disp(['connect the terminals ' num2str(rows(i)) ' and ' num2str(cols(i)) ])
    end
end

% MODIFICATION LOG
%
% 051122 med   Created
% 060116 per   Added documentation
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
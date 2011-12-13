%% Network Reliability
%
%% Problem description
% We consider the military telecommunications network represented
% below. This network consists of eleven sites connected by
% bidirectional lines for data transmission. For reliability reasons
% in the case of a conflict, the specifications require that the two
% sites (nodes) 10 and 11 of the network remain able to communicate
% even if any three other sites of the network are destroyed. Does
% the network satisfy this requirement?
%
%
%         1 --- 2 ---------- 8
%
%      /  |   / |          / |
%     /   |  /  |         /  |
%    /    | /   |        /   |
%               |            |
%  11 --- 3  ---+---- 10     |
%               |            |
%   | \ /   \   |   / |  \   |
%   |  X     \  |  /  |   |  |
%   | / \     \ | /   |   |  |
%                     |   |  |
%   4 --- 5 --- 9     |   |  |
%                     |   |  |
%     \         |    /    |  |
%      \        |   /     \  |
%       \       |  /       \ |
%        \
%         ----- 6 ---------- 7
%
%
%
%% Variables
%
%  arcs_out/in        For an arc i arcs_out(i) is its starting node
%                     and arcs_in(i) ts target node
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
arcs_in    = [1 1  1 2 2 2 3 3  3  3 4 4  4 5  5 6 6  6 7  7  8  9]';
arcs_out   = [2 3 11 3 8 9 4 9 10 11 5 6 11 9 11 7 9 10 8 10 10 10]';

n1 = length(arcs_in);   % arcs_in
n  = 2*n1;              % bi-directional
nodes = length(unique([arcs_in; arcs_out]));
arcs_in_new = [arcs_in;arcs_out];
arcs_out = [arcs_out;arcs_in];
arcs_in  = arcs_in_new;

arcs = tom('arcs',n,1,'int');

% All variables are binary.
bnds = {0 <= arcs <= 1};

% Kirchhoff's law, 1 source and 1 sink
con1 = cell(nodes-2,1);
for i=1:nodes-2 %Assuming 10 and 11 are source and sink
    idx1 = find(arcs_in == i);
    idx3 = find(arcs_out == i);
    con1{i} = {sum(arcs(idx1)) == sum(arcs(idx3))};
end

% Limit flow in intermediate node
con2 = cell(nodes-2,1);
for i=1:nodes-2 %Assuming 10 and 11 are source and sink
    idx1 = find(arcs_in == i);
    con2{i} = {sum(arcs(idx1)) <= 1};
end

% Source has no incoming flows
idx1 = find(arcs_out == nodes-1);
con3 = {sum(arcs(idx1)) == 0};

% Objective
idx1 = find(arcs_in == (nodes-1));
objective = -sum(arcs(idx1));

constraints = {bnds, con1, con2, con3};
options = struct;
options.solver = 'cplex';
options.name   = 'Network Reliability';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    f      = -subs(objective,sol);
    x      = -sol.arcs;
    disp(['The network can manage the loss of ' num2str(f-1) ' nodes.']);
    disp(['There are ' num2str(f) ' disjunctive paths.']);
    [arcs,direction] = find(reshape(x,length(arcs_in)/2,2));
    disp('For a maximal number of disjunctive paths, activate...');
    for arc = 1:length(arcs),
        if direction(arc) == 1,
            disp([' arc ' num2str([arcs_out(arcs(arc))]) '->' ...
                num2str([arcs_in(arcs(arc))])]);
        end
        if direction(arc) == 2,
            disp([' arc ' num2str([arcs_in(arcs(arc))]) '->'  ...
                num2str([arcs_out(arcs(arc))])]);
        end
    end
end

% MODIFICATION LOG
%
% 051107 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
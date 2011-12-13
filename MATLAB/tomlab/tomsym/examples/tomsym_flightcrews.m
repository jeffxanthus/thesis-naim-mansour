%% Composing Flight Crews
%
%% Problem description
% During the Second World War the Royal Air Force (RAF) had many
% foreign pilots speaking different languages and more or less
% trained on the different types of aircraft. The RAF had to form
% pilot/co-pilot pairs (‘crews’) for every plane with a compatible
% language and a sufficiently good knowledge of the aircraft type. In
% our example there are eight pilots. In the following table every
% pilot is characterized by a mark between 0 (worst) and 20 (best)
% for his language skills (in English, French, Dutch, and Norwegian),
% and for his experience with different two-seater aircraft
% (reconnaissance, transport, bomber, fighterbomber, and supply
% planes).
%
% Ratings of pilots
%
%  +----------+-------------------------+--+--+--+--+--+--+--+--+
%  |          |Pilot                    | 1| 2| 3| 4| 5| 6| 7| 8|
%  +----------+-------------------------+--+--+--+--+--+--+--+--+
%  |Language  |English                  |20|14| 0|13| 0| 0| 8| 8|
%  |          |French                   |12| 0| 0|10|15|20| 8| 9|
%  |          |Dutch                    | 0|20|12| 0| 8|11|14|12|
%  |          |Norwegian                | 0| 0| 0| 0|17| 0| 0|16|
%  +----------+-------------------------+--+--+--+--+--+--+--+--+
%  |Plane type|Reconnaissance           |18|12|15| 0| 0| 0| 8| 0|
%  |          |Transport                |10| 0| 9|14|15| 8|12|13|
%  |          |Bomber                   | 0|17| 0|11|13|10| 0| 0|
%  |          |Fighter-bomber           | 0| 0|14| 0| 0|12|16| 0|
%  |          |Supply plane             | 0| 0| 0| 0|12|18| 0|18|
%  +----------+-------------------------+--+--+--+--+--+--+--+--+
%
% A valid flight crew consists of two pilots that both have each at
% least 10/20 for the same language and 10/20 on the same aircraft
% type.
%
% Question 1:
% Is it possible to have all pilots fly?
%
% Subsequently, we calculate for every valid flight crew the sum of
% their scores for every aircraft type for which both pilots are
% rated at least 10/20. This allows us to define for every crew the
% maximum score among these marks. For example, pilots 5 and 6 have
% marks 13 and 10 on bombers and 12 and 18 on supply planes. The
% score for this crew is therefore max(13 + 10, 12 + 18) = 30.
%
% Question 2:
% Which is the set of crews with maximum total score?
%
% Compatibility graph for pilots
% Nodes = pilots
% Arcs  = possible combinations (with scores)
%
%   1 ----30--- 2 ----27--- 3
%
%   | \       / |         / |
%   |  24   28  |        /  |
%   |   \   /   27     26   |
%   |           |      /    |
%   |     4     |     /    30
%  25           |    /      |
%   |   /   \   |   /       |
%   |  29    21 |  /        |
%   | /       \ | /         |
%
%   5 ---30---- 6 -----28-- 7
%
%     \         |         /
%      \        |        /
%       \       |       /
%        \      36     /
%         \     |     /
%         30    |   25
%           \   |   /
%            \  |  /
%             \ | /
%
%               8
%
%
%% Variables
%
%  scores                     Scores for language and piloting skills
%  arcs_out/in                For an arc i arcs_out(i) is its starting node
%                             and arcs_in(i) ts target node
%  arcs_score                 Strength of an arc
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
scores = [ 20 14  0 13  0  0  8  8;...
    12  0  0 10 15 20  8  9;...
    0 20 12  0  8 11 14 12;...
    0  0  0  0 17  0  0 16;...
    18 12 15  0  0  0  8  0;...
    10  0  9 14 15  8 12 13;...
    0 17  0 11 13 10  0  0;...
    0  0 14  0  0 12 16  0;...
    0  0  0  0 12 18  0 18];

arcs_in    = [1 1 1 2 2 2 3 3 4 4 5 5 6 6 7]';
arcs_out   = [5 4 2 4 6 3 6 7 5 6 6 8 8 7 8]';

arcs_score = [25 24 30 28 27 27 26 30 29 21 30 30 ...
    36 28 25]';

n = length(arcs_in);
n1 = size(scores,2);

fly = tom('fly',n,1,'int');

% All variables are integers.
bnds = {0 <= fly <= 1};

% Orig constr.
con = cell(n1,1);
for i=1:n1
    idx_in  = find(arcs_in  == i);
    idx_out = find(arcs_out == i);
    idx = unique([idx_in;idx_out]);
    con{i} = {sum(fly(idx)) <= 1};
end

% Objective
objective = -arcs_score'*fly;

constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Composing Flight Crews';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    idx = find(sol.fly);
    for crew = 1:idx,
        id = idx(crew);
        disp(['crew ' num2str(crew) ' has pilots ' ...
            num2str(arcs_in(id)) ' and ' num2str(arcs_out(id))])
    end
end

% MODIFICATION LOG
%
% 051021 med   Created
% 060113 per   Added documentation
% 060125 per   Moved disp to end
% 060203 med   Updated print level
% 090308 med   Converted to tomSym
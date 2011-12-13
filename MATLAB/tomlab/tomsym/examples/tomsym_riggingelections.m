%% Rigging Elections
%
%% Problem description
% In a country not very far away, the party of the duke Sark Mevo has
% finally beaten the ruling coalition headed by the princess Reguel
% Tekris. Mevo wishes to consolidate his position in the capital, the
% fourteen quarters of which need to be grouped to electoral
% districts. A schematic map of the capital is given in the figure
% below. The quarters are numbered from 1 to 14. The two other
% numbers are the forecast number of favorable votes for Mevo and the
% total number of electors per quarter. All electors must vote and
% the winner needs to have the absolute majority. A valid electoral
% district is formed by several adjacent quarters and must have
% between 30,000 and 100,000 voters. Two quarters that touch each
% other just at a point like 12 and 13 are not considered adjacent.
% Electoral districts consisting of a single quarter are permitted if
% it has at least 50,000 voters. Nevertheless, Mevo may not decently
% define an electoral district solely with quarter 10, since this
% contains his residence. Determine a partitioning into five
% electoral districts that maximizes the number of seats for Mevo.
% Should this cause any difficulties, try a partitioning into six
% districts. Snirp, the mathematical jester, suggests Mevo uses
% Mathematical Programming...
%
%  +-------------+--+------+-------------------------------+
%  |    1        |  |   6  |  7  12000/30000               |
%  | 17500/30000 |  | 9000/+---------------+---------------+
%  +-------------+  |40000 |    8          |               |
%  |    2        |  |      |  10000/30000  |  9            |
%  | 15000/50000 |  +------+-------+-------+ 26000/40000   |
%  +-------------+    5    |       |  11   |               |
%  |    3        |  18000/ |  10   | 2500/ +---------------+
%  | 14200/20000 |  20000  | 34000/| 10000 | 12 27000/60000|
%  +-------------+---------+ 60000 +-------+---------------+
%  |         4             |       |  13   |  14           |
%  |      42000/70000      |       | 29000/| 15000/40000   |
%  +                       |       | 40000 |               |
%  +-----------------------+-------+-------+---------------+
%
% Map of the capital and its quarters.
% Legend: quarter number, supporters/electorate
%
%% Variables
%
%  in/out              Quarter in(i) borders to quarter out(i)
%  population          Population of each quarter
%  votes               Supporters er quarter
%  minpop              Min size of an electoral district
%  maxpop              Max size of an electoral district
%  minsingle           Min size of a quarter to be a district
%  districtsnum        Number of districts wanted
%  illegalsingle       This quarter may not be single
%  districts           The possible partition of the quarters
%                      into desired number of districts.
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
function tomsym_riggingelections

in  =  [1 1 2 2 3 3 4  4 5  5 6 6 7 7 8  8  8  9  9 10 10 11 11 12 13]';
out =  [2 5 3 5 4 5 5 10 6 10 7 8 8 9 9 10 11 11 12 11 13 12 13 14 14]';

population = [30 50 20 70 20 40 30 30 40 60 10 60 40 40]';
votes      = [17.5 15 14.2 42 18 9 12 10 26 34 2.5 27 29 15]';
minpop     = 30;
maxpop     = 100;
minsingle  = 50;
districtsnum = 6;
illegalsingle = 10;

districts = [];
rowlen = length(population);
maxq   = rowlen;

for q=1:maxq
    if (population(q) >= minsingle & q ~= illegalsingle)
        districts = [districts; zeros(1,maxq)];
        districts(end, q) = 1;
    end
    idx = find(in == q);
    for i = 1:length(idx)
        p = out(idx(i));
        if (population(q) + population(p) <= maxpop)
            if (population(q) + population(p) >= minpop)
                districts = addNeighbor(districts, q, p, rowlen);
                idx2 = find(in == p);
                for j = 1:length(idx2)
                    r = out(idx2(j));
                    if (population(q) + population(p) + population(r) <= maxpop)
                        if (population(q) + population(p) + population(r) >= minpop)
                            districts = addNeighbor(districts, q, p, rowlen, r);
                        end
                    end
                end
            end
        end
    end
end

n = size(districts, 1);  %possible districts

% Calculate majority
d = zeros(n,1);
for i=1:n
    idx = find(districts(i,:) == 1);
    if sum(votes(idx))/sum(population(idx)) > 0.5
        d(i,1) = 1;
    end
end

choose = tom('choose',n,1,'int');

% All variables are binary
bnds = {0 <= choose <= 1};

% Districts chosen
con1 = {sum(choose) == districtsnum};

% Quarter only once
con2 = {districts'*choose == 1};

% Objective
objective = -d'*choose;
constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Rigging Elections';
sol = ezsolve(objective,constraints,[],options);

sol.districts = districts;

PriLev = 1;
if PriLev > 0
    temp   = sol.districts(find(sol.choose),:);
    disp('to rig the elections:')
    for d  = 1:size(temp,1),
        qs = find(temp(d,:));
        disp(['   merge the quarters ' num2str(qs) ' into a district.'])
    end
end

function districts = addNeighbor(districts, q, p, rowlen, r)
if nargin <5
    districts = [districts; zeros(1,rowlen)];
    districts(end, [q, p]) = [1 1];
else
    districts = [districts; zeros(1,rowlen)];
    districts(end, [q, p, r]) = [1 1 1];
end

% MODIFICATION LOG
%
% 051205 med   Created
% 060118 per   Added documentation
% 060125 per   Moved disp to end
% 090325 med   Converted to tomSym
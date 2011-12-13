%% Alcuin's River Crossing
% TomSym implementation of GAMS Example (CROSS,SEQ=191)
%
% A farmer carrying a bushel of corn and accompanied by a goose and a wolf
% came to a river. He found a boat capable of transporting himself plus one
% of his possessions - corn, goose, or wolf - but no more. Now, he couldn't
% leave the corn alone with the goose, nor the goose alone with the wolf,
% else one would consume the other. Nevertheless, he succeeded in getting
% himself and his goods across the river safely.
%
% Borndoerfer, R, Groetschel, M, and Loebel, A, Alcuin's
% Transportation Problem and Integer Programming. Konrad Zuse
% Zentrum for Informationstechnik, Berlin, 1995.
%
% Contributed by Soren Nielsen, Institute for Mathematical Sciences
%                University of Copenhagen
%
% i: items (goose, wolf, corn)
%
% t: time  (t1-t10);

% Crossing - near to far is +1 - far to near -1
dirt = (-1).^((1:10)'-1);

% 1 if the item is on the far side at time t
toms 3x10 integer y % binary

% crossing the river
rivercross = tom('rivercross',3,10); %toms 3x10 rivercross
toms 10x1 done % all items on far side
toms nocross   % number of non crossing periods

cbnd = {0 <= y <= 1 % y is binary, rest positive
    rivercross >= 0
    done >= 0
    nocross >= 0};

% Cross definition
eq1 = {};
for i=1:3
    for t=1:9
        eq1 = {eq1; y(i,t+1) == y(i,t) + dirt(t)*rivercross(i,t)};
    end
end

% Everything on far side
eq2 = {};
for i=1:3
    for t=1:10
        eq2 = {eq2; done(t) <= y(i,t)};
    end
end

% Limit rivercross
eq3 = {};
for t=1:9
    eq3 = {eq3; sum(rivercross(:,t)) <= 1};
end

% Eat none 1
eq4 = {};
for t=1:10
    eq4 = {eq4; dirt(t)*(y(1,t) + y(2,t) - 1) <= done(t)};
end

% Eat none 2
eq5 = {};
for t=1:10
    eq5 = {eq5; dirt(t)*(y(1,t) + y(3,t) - 1) <= done(t)};
end

% Objective
nocross = sum(done);

cbndinit = {y(:,1) == 0};

con1 = {cbnd; eq1; eq2; eq3; eq4; eq5; cbndinit};
solution = ezsolve(-nocross,con1);
%% Job Shop Scheduling
%
%% Problem description
% A company has received an order for three types of wallpapers: one
% (paper 1) has a blue background with yellow patterns, another
% (paper 2) a green background and blue and yellow patterns, and the
% last (paper 3) has a yellow background with blue and green
% patterns. Every paper type is produced as a continuous roll of
% paper that passes through several machines, each printing a
% different color. The order in which the papers are run through the
% machines depends on the design of the paper: for paper 1 first the
% blue background and then the yellow pattern is printed. After the
% green background for paper 2, first the blue and then the yellow
% patterns are printed. The printing of paper 3 starts with the
% yellow background, followed by the blue and then the green
% patterns.
%
%       p2             p1              p3
%       |              |               |
%       V              V               V
%      +-----+        +----+ --p1-> +------+
%      |GREEN| --p2-> |BLUE| --p2-> |YELLOW|
%      +-----+ <-p3-- +----+ <-p3-- +------+
%
%         |                           |   |
%         V                           V   V
%         p3                          p1  p2
%
% Production flows through printing machines
%
% The processing times differ depending on the surface that needs to
% be printed. The times (in minutes) for applying every color of the
% three paper types are given in the following table.
%
% Times required for applying every color
%
%  +-------+------+-------+-------+-------+
%  |Machine|Color |Paper 1|Paper 2|Paper 3|
%  +-------+------+-------+-------+-------+
%  |   1   |Blue  |  45   |  20   |  12   |
%  |   2   |Green |   -   |  10   |  17   |
%  |   3   |Yellow|  10   |  34   |  28   |
%  +-------+------+-------+-------+-------+
%
% Knowing that every machine can only process one wallpaper at a time
% and that a paper cannot be processed by several machines
% simultaneously, how should the paper printing be scheduled on the
% machines in order to finish the order as early as possible?
%
%% Variables
%
% the colors are ordered as follows:
%
%  color1 = blue   (paper 1's first color)
%  color2 = green  (paper 2's first color)
%  color3 = yellow (paper 3's first color)
%
%  proctimes         Times for the papers in the machines
%  flows             The order of colors on the paperS
%  finalm            The last machine for each paper
%  bigM              The total processingtimes.
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
proctimes = [45 20 12 10 17 10 34 28]';

flows = [1 3 0;...
    2 1 3;...
    3 1 2];

finalm = [3;3;2];

bigM  = sum(sum(proctimes));

disj1 = [1 1 2 4 6 6 7]';
disj2 = [2 3 3 5 7 8 8]';

arc1 = [1 4 2 8 3]';
arc2 = [6 2 7 3 5]';

tasks = 8;

start  = tom('start',tasks,1);
y      = tom('y',length(disj1),1,'int');
finish = tom('finish',1,1); %Scalar

% All slots are integers
bnds1 = {start >= 0, 0 <= finish <= bigM, 0 <= y <= 1};

% Task constraint
taskcon = {start(1:8) + proctimes(1:8) <= finish};

% Arc constraint
arccon = {};
for i=1:length(arc1)
    arccon{i} = {start(arc1(i)) + proctimes(arc1(i)) <= start(arc2(i))};
end

% Disj constraint
disjcon = {};
for i=1:length(disj1)
    disjcon{i} = ...
        {start(disj1(i)) + proctimes(disj1(i)) <= ...
        start(disj2(i)) + bigM*y(i)};
    disjcon{i+length(disj1)} = ...
        {start(disj2(i)) + proctimes(disj2(i)) <= ...
        start(disj1(i)) + bigM*(1-y(i))};
end

constraints = {bnds1, taskcon, arccon, disjcon};
options = struct;
options.solver = 'cplex';
options.name   = 'Job Shop Scheduling';
sol = ezsolve(finish,constraints,[],options);

PriLev = 1;
if PriLev > 0
    c      = 3; % number of colors
    p      = 3; % number of papers

    disp(['all papers are finished after ' num2str(sol.finish) ' minutes'])
    disp(' ')
    for j = 1:p,
        disp(['paper ' num2str(j) ':'])
        if j==1
            colortimes = sol.start(1:3);
        elseif j == 2
            colortimes = [0;sol.start(4:5)];
        else
            colortimes = sol.start(6:8);
        end

        [time,col] = sort(colortimes);
        for i = 1:c,
            this_col  = col(i);
            this_time = time(i);
            if this_time == 0
                disp(['   starts using color ' num2str(this_col) ' after ' ...
                    num2str(this_time) ' minutes (or not at all)'])
            else
                disp(['   starts using color ' num2str(this_col) ' after ' ...
                    num2str(this_time) ' minutes'])
            end
        end
    end
end

% MODIFICATION LOG
%
% 060102 med   Created
% 060111 per   Added documentation
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym
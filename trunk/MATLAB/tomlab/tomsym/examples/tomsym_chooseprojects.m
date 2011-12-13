%% Choice of Expansion Projects
%
%% Problem description
% The large company Tatayo in the north of Italy has specialized in
% the construction of cars for more than ten years. The company
% wishes to expand and has issued internally a call for proposals for
% expansion projects for a planning period of five years. Among the
% many, often cranky, propositions the management has retained five
% projects. Every project has an annual cost and is designed to
% produce a benefit after five years. The first table below gives a
% list of the projects with short descriptions and the expected
% benefit after five years. The forecast annual costs of the projects
% for the next five years are detailed in the second table below,
% together with the funds available. Which project(s) should the
% management choose now to maximize the total benefit after five
% years?
%
% Estimated benefits of the projects (in million $)
%
%  +-------+------------------------------+----------------+
%  |Project|Description                   |Expected benefit|
%  +-------+------------------------------+----------------+
%  |   1   |Expand assembly line          |     10.8       |
%  |   2   |Reorganize the main shop      |      4.8       |
%  |   3   |New painting facilities       |      3.2       |
%  |   4   |Research for a new concept car|      4.44      |
%  |   5   |Reorganize the logistics chain|     12.25      |
%  +-------+------------------------------+----------------+
%
% Annual costs of projects and available funds (in million $)
%
%  +-------+------+------+------+------+------+
%  |Project|Year 1|Year 2|Year 3|Year 4|Year 5|
%  +-------+------+------+------+------+------+
%  |   1   | 1.8  | 2.4  | 2.4  | 1.8  | 1.5  |
%  |   2   | 1.2  | 1.8  | 2.4  | 0.6  | 0.5  |
%  |   3   | 1.2  | 1.0  | 0.0  | 0.48 | 0.0  |
%  |   4   | 1.4  | 1.4  | 1.2  | 1.2  | 1.2  |
%  |   5   | 1.6  | 2.1  | 2.5  | 2.0  | 1.8  |
%  +-------+------+------+------+------+------+
%  |Funds  | 4.8  | 6.0  | 4.8  | 4.2  | 3.5  |
%  +-------+------+------+------+------+------+
%
%% Variables
%
%  benefit                    Expected benefit
%  budget                     Funds available each year
%  costmat                    Cost per project and year
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
benefit       = [10.8 4.8 3.2 4.44 12.25]'*1e6;
budget        = [ 4.8 6.0 4.8  4.2   3.5]'*1e6;

costmat       = [1.8 2.4 2.4 1.8 1.5;...
    1.2 1.8 2.4 0.6 0.5;...
    1.2 1.0 0.0 .48 0.0;...
    1.4 1.4 1.2 1.2 1.2;...
    1.6 2.1 2.5 2.0 1.8]*1e6;

n = length(benefit);   %projects

choose = tom('choose',n,1,'int');

% All variables are integer.
bnds = {0 <= choose <= 1};

% Cost constraints
con = {(choose'*costmat)' <= budget};

% Objective
objective = -benefit'*choose;

constraints = {bnds, con};
options = struct;
options.solver = 'cplex';
options.name   = 'Choice of Expansion Projects';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    names  = ['  Expand assembly line          ' ;
        '  Reorganize the main shop      ' ;
        '  New painting facilities       ' ;
        '  Research for a new concept car' ;
        '  Reorganize the logistics chain'];
    idx    = find(sol.choose);
    disp('The management should choose the following projects:')
    disp(names(idx,:))
end

% MODIFICATION LOG
%
% 051201 med   Created
% 060117 per   Added documentation
% 090308 med   Converted to tomSym
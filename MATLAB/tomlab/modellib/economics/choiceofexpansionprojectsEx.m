% function Result = choiceofexpansionprojectsEx(PriLev)
%
% Creates a TOMLAB MIP problem for choice of expansion projects
%
% CHOICE OF EXPANSION PROJECTS
%
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
% +-------+------------------------------+----------------+
% |Project|Description                   |Expected benefit|
% +-------+------------------------------+----------------+
% |   1   |Expand assembly line          |     10.8       |
% |   2   |Reorganize the main shop      |      4.8       |
% |   3   |New painting facilities       |      3.2       |
% |   4   |Research for a new concept car|      4.44      |
% |   5   |Reorganize the logistics chain|     12.25      |
% +-------+------------------------------+----------------+
% 
% Annual costs of projects and available funds (in million $)
% 
% +-------+------+------+------+------+------+
% |Project|Year 1|Year 2|Year 3|Year 4|Year 5|
% +-------+------+------+------+------+------+
% |   1   | 1.8  | 2.4  | 2.4  | 1.8  | 1.5  |
% |   2   | 1.2  | 1.8  | 2.4  | 0.6  | 0.5  |
% |   3   | 1.2  | 1.0  | 0.0  | 0.48 | 0.0  |
% |   4   | 1.4  | 1.4  | 1.2  | 1.2  | 1.2  |
% |   5   | 1.6  | 2.1  | 2.5  | 2.0  | 1.8  |
% +-------+------+------+------+------+------+
% |Funds  | 4.8  | 6.0  | 4.8  | 4.2  | 3.5  |
% +-------+------+------+------+------+------+
%
% VARIABLES
% 
% benefit                    Expected benefit
% budget                     Funds available each year
% costmat                    Cost per project and year
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = choiceofexpansionprojectsEx(2);
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
% Written Dec 1, 2005.   Last modified Dec 1, 2005.

function Result = choiceofexpansionprojectsEx(PriLev)

if nargin < 1
   PriLev = 1;
end

benefit       = [10.8 4.8 3.2 4.44 12.25]'*1e6;
budget        = [ 4.8 6.0 4.8  4.2   3.5]'*1e6;

costmat       = [1.8 2.4 2.4 1.8 1.5;...
      1.2 1.8 2.4 0.6 0.5;...
      1.2 1.0 0.0 .48 0.0;...
      1.4 1.4 1.2 1.2 1.2;...
      1.6 2.1 2.5 2.0 1.8]*1e6;

Prob = choiceofexpansionprojects(benefit, budget, costmat);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   names  = ['  Expand assembly line          ' ;
      '  Reorganize the main shop      ' ;
      '  New painting facilities       ' ;
      '  Research for a new concept car' ;
      '  Reorganize the logistics chain'];
   idx    = find(Result.x_k);
   disp('The management should choose the following projects:')
   disp(names(idx,:))
end

% MODIFICATION LOG
%
% 051201 med   Created.
% 060117 per   Added documentation.

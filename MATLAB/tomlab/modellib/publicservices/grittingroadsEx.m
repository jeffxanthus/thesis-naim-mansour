% function Result = grittingroadsEx(PriLev)
%
% Creates a TOMLAB MILP problem for gritting roads
%
% GRITTING ROADS
%
% In the case of ice, all the streets of a village need to be
% gritted. A schematic map of the streets is given in the figure
% below. The nodes correspond to street intersections, the arcs to
% the streets that need to be gritted. The arcs are labeled with the
% length of the street in meters.
%     
% Graph of the village streets
%
%      --150->   --130->   --100->
%    1         2         3         4
%                <-140--   <-100--
%   | ^      /| ^        ^        | ^
%   | 165   / | 170      |        | 180
%   | |    /  | |       200       | |
%  165|  230 160|        |       190|
%   | |  /    | |        |        | |
%   | | /     | |        |        | |
%   V |V      V |        |        V |
%      --144->   --128->      
%    5         6         7 --109-> 8 
%      <-144--   <-122--      
%    ^       /| ^      ^| ^      / | 
%    194    / |174    / | |     /  | 
%    |     /  | |    / 185|185 /   | 
%    |   218 174|  233  | |   /   190
%    |   /    | |  /    | |  140   | 
%    |  /     | | /     | | /      | 
%    | V      V |/      V |V       V 
%                
%    9 --148-> 10 <-135- 11 -110-> 12 
%                      
% The highway maintenance depot with the gritting truck is located at
% intersection number 1. The truck has a sufficiently large capacity
% to grit all the streets during a single tour. Due to the one-way
% streets, it may be forced to pass several times through the same
% street that has already been gritted. Determine a tour for the
% gritting truck of minimum total length that passes through all
% streets. For bidirectional streets, both directions need to be
% gritted separately.
%
% VARIABLES
%
% in/out                     Road i starts in in(i) and 
%                            goes out to out(i)
% lengths                    Lengths of the roads
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = grittingroadsEx(2);
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
% Written Dec 6, 2005.   Last modified Dec 6, 2005.

function Result = grittingroadsEx(PriLev)

if nargin < 1
   PriLev = 1;
end

in  =  [1 1 2 2 2 3 3 4 4 5 5 6 6 6 6  6 7 7 7  7 8  8  8 9 9 10 10 11 11 12]';
out =  [2 5 3 5 6 2 4 3 8 1 6 2 5 7 9 10 3 6 8 11 4 11 12 5 10 6  7  7 10 11]';

lengths = [150 165 130 230 160 140 100 100 190 165 144 170 144 128 218 174 ...
      200 122 109 185 180 141 190 194 148 174 233 185 135 110]';

Prob = grittingroads(in, out, lengths);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   twice  = find(Result.x_k > 1);
   disp('all roads travelled once, except:')
   for i  = 1:length(twice),
      idx  = twice(i);
      disp(['  road from ' num2str(in(idx)) ' to ' num2str(out(idx)) ...
            ' ( that is travelled ' num2str(Result.x_k(idx))   ' times)'])
   end
end

% MODIFICATION LOG
%
% 051206 med   Created.
% 060118 per   Added documentation.
% 060125 per   Moved disp to end

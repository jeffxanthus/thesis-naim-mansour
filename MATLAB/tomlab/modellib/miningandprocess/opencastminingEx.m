% function Result = opencastminingEx(PriLev)
%
% Creates a TOMLAB MIP problem for open cast mining
%
% OPENCAST MINING
%
% PROBLEM
%
% An opencast uranium mine is being prospected. Based on the results of
% some test drillings the mine has been subdivided into exploitation units
% called blocks. The pit needs to be terraced to allow the trucks to drive
% down to its bottom. The uranium deposit extends from east to west. The 
% pit is limited in the west by a village and in the east by a group of
% mountains. Taking into account these constraints, 18 blocks of 10,000 
% tonnes on three levels have been identified (Figure 6.3). To extract a
% block, three blocks of the level above it need to be extracted: the
% block immediately on top of it, and also, due to the constraints on the
% slope, the blocks to the right and to the left.
%
%
%   Village                                 Mountains
%          +---+---+---+---+---+---+---+---+
% Level 1: | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
%          +---+---+---+---+---+---+---+---+
% Level 2: |   | 9 | 10| 11| 12| 13| 14|   |
%          +---+---+---+---+---+---+---+---+
% Level 3: |   |   | 15| 16| 17| 18|   |   |
%          +---+---+---+---+---+---+---+---+
%
% For example: If we were to extract block 17 we also need to extract the
%              blocks 11, 12 and 13. And in order to extract these blocks
%              we also need to extract blocks 3, 4, 5, 6 and 7.
%
% It costs $ 100 per tonne to extract a block of level 1, $ 200 per
% tonne for a block of level 2, and $ 300 per tonne for a block of
% level 3, with the exception of the hatched blocks that are formed of a
% very hard rock rich in quartz and cost $ 1000 per ton. The only blocks
% that contain uranium are those displayed in a gray shade (1, 7, 10, 12,
% 17, 18). Their market value is 200, 300, 500, 200, 1000, and 
% $ 1200/tonne respectively. Block 18, although rich in ore, is made of 
% the same hard rock as the other hatched blocks. Which blocks should be 
% extracted to maximize the total benefit?
%
% VARIABLES
%
% values                     The profit for extracting the blocks
% depends                    Block dependencies. 
%
% RESULT
%
% The result will be a vector with binary variables, as many as we have 
% blocks. A 1 means do extract, and a 0 means do not extract.
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
% Result       Result structure.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 7, 2005.   Last modified Oct 7, 2005.

function Result = opencastminingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

values      = [200-100;-100;-100;-100;-100;-100;300-100;-100;...
      -1000;500-200;-200;200-200;-200;-1000;-1000;-1000;...
      1000-300;1200-1000]*10000;
depends     = [-1 -1 -1  0  0  0  0  0  3  0  0  0  0  0  0  0  0  0;...
      0 -1 -1 -1  0  0  0  0  0  3  0  0  0  0  0  0  0  0;...
      0  0 -1 -1 -1  0  0  0  0  0  3  0  0  0  0  0  0  0;... 
      0  0  0 -1 -1 -1  0  0  0  0  0  3  0  0  0  0  0  0;... 
      0  0  0  0 -1 -1 -1  0  0  0  0  0  3  0  0  0  0  0;... 
      0  0  0  0  0 -1 -1 -1  0  0  0  0  0  3  0  0  0  0;... 
      0  0  0  0  0  0  0  0 -1 -1 -1  0  0  0  3  0  0  0;... 
      0  0  0  0  0  0  0  0  0 -1 -1 -1  0  0  0  3  0  0;... 
      0  0  0  0  0  0  0  0  0  0 -1 -1 -1  0  0  0  3  0;... 
      0  0  0  0  0  0  0  0  0  0  0 -1 -1 -1  0  0  0  3];

Prob = opencastmining(values, depends);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   level1 = [];
   for i = 1 : 8,
      if Result.x_k(i) == 1,
         level1 = [level1 ' * '];
      else
         level1 = [level1 ' - '];
      end
   end
   level2 = [ ' - '];
   for i = 9 : 14,
      if Result.x_k(i) == 1,
         level2 = [level2 ' * '];
      else
         level2 = [level2 ' - '];
      end
   end
   level2 = [level2 ' - '];
   level3 = [ ' - ' ' - '];
   for i = 15 : 18,
      if Result.x_k(i) == 1,
         level3 = [level3 ' * '];
      else
         level3 = [level3 ' - '];
      end
   end
   level3 = [level3 ' - ' ' - '];
disp(' ')
disp('Extraction profile')
disp(' ')
disp(['Level 1:  ' level1])   
disp(['Level 2:  ' level2])   
disp(['Level 3:  ' level3])   
disp(' ')
disp('Legend: * means do extract')
disp('        - means do not extract')
end

% MODIFICATION LOG
%
% 051007 med   Created.
% 060109 per   Added documentation.
% 060125 per   Moved disp to end

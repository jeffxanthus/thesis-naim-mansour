% function Result = jobshopschedulingEx(PriLev)
%
% Creates a TOMLAB MIP problem for job shop scheduling
%
% JOB SHOP SCHEDULING
%
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
% +-------+------+-------+-------+-------+
% |Machine|Color |Paper 1|Paper 2|Paper 3|
% +-------+------+-------+-------+-------+
% |   1   |Blue  |  45   |  20   |  12   |
% |   2   |Green |   -   |  10   |  17   |
% |   3   |Yellow|  10   |  34   |  28   |
% +-------+------+-------+-------+-------+
%
% Knowing that every machine can only process one wallpaper at a time
% and that a paper cannot be processed by several machines
% simultaneously, how should the paper printing be scheduled on the
% machines in order to finish the order as early as possible?
%
% VARIABLES
%
% the colors are ordered as follows: 
%    color1 = blue   (paper 1's first color)
%    color2 = green  (paper 2's first color)
%    color3 = yellow (paper 3's first color)
%    
%
% proctimes                  Times for the papers in the machines
% flow                       The order of colors on the paperS
% final                      The last machine for each paper
% bigM                       The total processingtimes.
%
% RESULTS
%
% run this for explanation of x_k
% Result = jobshopschedulingEx(2);
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
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Jan 2, 2006.   Last modified Jan 2, 2006.

function Result = jobshopschedulingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

proctimes   = [45 20 12;...
      0  10 17;...
      10 34 28];

flow        = [1 3 0;...
      2 1 3;...
      3 1 2];

final       = [3;3;2];            

bigM        = sum(sum(proctimes));

Prob = jobshopscheduling(proctimes, flow, final, bigM);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   c      = 3; % number of colors
   p      = 3; % number of papers
   temp   = reshape(Result.x_k,c,c+c+2);
   disp(['all papers are finished after ' num2str(max(temp(:,c+1))) ' minutes'])
   disp(' ')
   for j = 1:p,
      disp(['paper ' num2str(j) ':'])
      colortimes = temp(j,1:c);
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
% 060102 med   Created.
% 060111 per   Added documentation.
% 060126 per   Moved disp to end

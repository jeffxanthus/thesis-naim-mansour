% function Result = efficiencyofhospitalsEx(PriLev)
%
% Creates a TOMLAB MILP problem for efficiency of hospitals
%
% EFFICIENCY OF HOSPITALS
%
% The administration of the hospitals in Paris decides to measure the
% efficiency of the surgery departments in four major hospitals with
% a desire to improve the service to the public. To keep this study
% anonymous, the hospitals are named H1 to H4. The method suggested
% to measure the efficiency is DEA (Data Envelopment Analysis). This
% method compares the performance of a fictitious hospital with the
% performances of the four hospitals.
%
% Three initial indicators (resources) are taken into account: the
% number of non-medical personnel, the general expenses, and the
% available number of beds. In addition, four final indicators
% (services) are analyzed: the number of hospital admissions per day,
% the number of consultations in the outpatients’ clinic, the number
% of nurses on duty every day, and the number of interns and doctors
% on duty every day. The corresponding data have been analyzed over a
% period of two years and the numbers representing a day of average
% activity in every hospital are given in the following two tables.
%
% Resource indicators
% 
% +---------------------+-----+------+-----+-----+
% |                     | H1  | H2   | H3  | H4  |
% +---------------------+-----+------+-----+-----+
% |Non-medical personnel| 90  | 87   | 51  | 66  |
% |General expenses (k$)|38.89|109.48|40.43|48.41|
% |Number of beds       | 34  | 33   | 20  | 33  |
% +---------------------+-----+------+-----+-----+
%
% Service indicators
%
% +-------------------+-----+-----+-----+-----+
% |                   | H1  | H2  | H3  | H4  |
% +-------------------+-----+-----+-----+-----+
% |Admissions         |30.12|18.54|20.88|10.42|
% |Consultations      |13.54|14.45| 8.52|17.74|
% |Interns and doctors| 13  |  7  |   8 | 26  |
% |Nurses on duty     | 79  | 55  |  47 | 50  |
% +-------------------+-----+-----+-----+-----+
% 
% Justify through the DEA method how hospital H2 is performing
% compared to the others.
%
% VARIABLES
%
% resources                  A matrix describing the resources
% services                   A matrix describing the services
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = efficiencyofhospitalsEx(2);
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

function Result = efficiencyofhospitalsEx(PriLev)

if nargin < 1
   PriLev = 1;
end

resources  = [90 87 51 66;...
      38.89 109.48 40.43 48.41;...
      34 33 20 33];
services   = [30.12 18.54 20.88 10.42;...
      13.54 14.45 8.52 17.74;...
      13 7 8 26;...
      79 55 47 50];

indices    = zeros(size(resources,2),1);

for i=1:size(resources,2)           
   Prob = efficiencyofhospitals(resources, services, i);
   Result = tomRun('cplex', Prob, PriLev);
   indices(i,1) = Result.x_k(end,1);
   Result.indices = indices;
end

if PriLev > 1,
   temp   = Result.indices;
   for i = 1:length(temp),
      disp(['H' num2str(i) ' is ' num2str(100*temp(i)) '% efficient'])
   end 
end


% MODIFICATION LOG
%
% 051206 med   Created.
% 060118 per   Added documentation.
% 060125 per   Moved disp to end

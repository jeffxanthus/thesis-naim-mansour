% function Result = examschedulingEx(PriLev)
%
% Creates a TOMLAB MILP problem for exam scheduling
%
% EXAM SCHEDULING
%
% At a technical university every term the third-year students choose
% eight modules from the eleven modules that are taught, depending on
% the option they wish to choose in the fourth year (there are two
% possible choices: "Production planning" and "Quality and security
% management"). In the current term, certain modules are obligatory
% for students who wish to continue with one of these options. The
% obligatory courses are Statistics (S), Graph models and algorithms
% (GMA), Production management (PM), and Discrete systems and events
% (DSE). The optional modules are: Data analysis (DA), Numerical
% analysis (NA), Mathematical programming (MP), C++, Java (J), Logic
% programming (LP), and Software engineering (SE). 
%
% Incompatibilities between different exams
% 
% +---+---+---+---+---+---+---+---+---+---+---+---+
% |   | DA| NA|C++| SE| PM| J |GMA| LP| MP| S |DSE|
% +---+---+---+---+---+---+---+---+---+---+---+---+
% |DA | - | X | - | - | X | - | X | - | - | X | X |
% |NA | X | - | - | - | X | - | X | - | - | X | X |
% |C++| - | - | - | X | X | X | X | - | X | X | X |
% |SE | - | - | X | - | X | X | X | - | - | X | X |
% |PM | X | X | X | X | - | X | X | X | X | X | X |
% |J  | - | - | X | X | X | - | X | - | X | X | X |
% |GMA| X | X | X | X | X | X | - | X | X | X | X |
% |LP | - | - | - | - | X | - | X | - | - | X | X |
% |MP | - | - | X | - | X | X | X | - | - | X | X |
% |S  | X | X | X | X | X | X | X | X | X | - | X |
% |DSE| X | X | X | X | X | X | X | X | X | X | - |
% +---+---+---+---+---+---+---+---+---+---+---+---+
% 
% Mrs Edeetee needs to schedule the exams at the end of the term.
% Every exam lasts two hours. Two days have been reserved for the
% exams with the following time slices: 8:00-10:00, 10:15-12:15,
% 14:00-16:00, and 16:15-18:15. For every exam she knows the set of
% incompatible exams that may not take place at the same time because
% they have to be taken by the same students. These incompatibilities
% are summarized in the table above.
%
% Help Mrs Edeetee construct a timetable so that no student has more
% than one exam at a time.
%
% VARIABLES
%
% incompatmat                Matrix of incompatabilities
% slots                      Timeslots (2 days * 4 per day)
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = examschedulingEx(2);

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
% Written Dec 5, 2005.   Last modified Dec 5, 2005.

function Result = examschedulingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

incompatmat     = [ 0 1 0 0 1 0 1 0 0 1 1;...
      1 0 0 0 1 0 1 0 0 1 1;...
      0 0 0 1 1 1 1 0 1 1 1;...
      0 0 1 0 1 1 1 0 0 1 1;...
      1 1 1 1 0 1 1 1 1 1 1;...
      0 0 1 1 1 0 1 0 1 1 1;...
      1 1 1 1 1 1 0 1 1 1 1;...
      0 0 0 0 1 0 1 0 0 1 1;...
      0 0 1 0 1 1 1 0 0 1 1;...
      1 1 1 1 1 1 1 1 1 0 1;...
      1 1 1 1 1 1 1 1 1 1 0];

slots           = 8;             

Prob = examscheduling(incompatmat, slots);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   e_names = ['DA  '; 'NA  '; 'C++ '; 'SE  ';
      'PM  '; 'J   '; 'GMA '; 'LP  ';  
      'MP  '; 'S   '; 'DSE '] ;
   exams   = length(e_names); % number of exams
   days    = 2;               % number of days
   times   = ['0800-1000'; '1015-1215'; '1400-1600'; '1615-1815'; ];
   slots   = size(times,1);   % time slots per day
   temp    = reshape(Result.x_k,exams,days*slots);
   
   
   for d = 1:days,
      disp([' day number ' num2str(d)])
      for s = 1:slots,
         i = s + (slots)*(d-1);
         exams_now = find(temp(:,i));
         e_list = [];
         for e = 1:length(exams_now),
            e_list = [e_list e_names(exams_now(e),:)];
         end
         if ~isempty(e_list),
            disp(['  during ' times(s,:) ': ' num2str(e_list)])
         end
      end
   end
end

% MODIFICATION LOG
%
% 051205 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end

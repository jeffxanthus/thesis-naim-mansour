% function Result = establishingacollegetimetableEx(PriLev)
%
% Creates a TOMLAB MILP problem for establishing a college time table
%
% ESTABLISHING A COLLEGE TIMETABLE
%
% Mr. Miller is in charge of establishing the weekly timetable for
% two classes of the last year in a college. The two classes have
% the same teachers, except for mathematics and sport. In the
% college all lessons have a duration of two hours. Furthermore,
% all students of the same class attend exactly the same courses.
% From Monday to Friday, the slots for courses are the following:
% 8:00–10:00, 10:15–12:15, 14:00–16:00, and 16:15–18:15. The
% following table lists the number of two-hour lessons that every
% teacher has to teach the students of the two classes per week.
%
% Number of 2-hour lessons per teacher and class
%+------------+-----------------+-------------------+-------------------+
%|Teacher     |Subject          |Lessons for class 1|Lessons for class 2|
%+------------+-----------------+-------------------+-------------------+
%|Mr Cheese   |English          |         1         |         1         |
%|Mrs Insulin |Biology          |         3         |         3         |
%|Mr Map      |History-Geography|         2         |         2         |
%|Mr Effofecks|Mathematics      |         0         |         4         |
%|Mrs Derivate|Mathematics      |         4         |         0         |
%|Mrs Electron|Physics          |         3         |         3         |
%|Mr Wise     |Philosophy       |         1         |         1         |
%|Mr Muscle   |Sport            |         1         |         0         |
%|Mrs Biceps  |Sport            |         0         |         1         |
%+------------+-----------------+-------------------+-------------------+
%
% The sport lessons have to take place on Thursday afternoon from
% 14:00 to 16:00. Furthermore, the first time slot on Monday morning
% is reserved for supervised homework. Mr Effofecks is absent every
% Monday morning because he teaches some courses at another college.
% Mrs Insulin does not work on Wednesday. And finally, to prevent
% students from getting bored, every class may only have one two-hour
% lesson per subject on a single day. Write a mathematical program
% that allows Mr Miller to determine the weekly timetable for the two
% classes.
%
% VARIABLES
%
% lessons1/2                 Number of lessons per subject and class
% subject                    Subject indices
% slots                      Possible slots
%
% RESULTS
%
% For an interpretation of the results, run:
% Result = establishingacollegetimetableEx(2);
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

function Result = establishingacollegetimetableEx(PriLev)

if nargin < 1
   PriLev = 1;
end

lessons1        = [ 1  3  2  0  4  3  1  1  0]';
lessons2        = [ 1  3  2  4  0  3  1  0  1]';
subject         = [ 1  2  3  4  4  5  6  7  7]';
slots           = 4*5;

Prob = establishingacollegetimetable(lessons1, lessons2, subject, slots);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   t_names = ['Cheese   '; 'Insulin  '; 'Map      ';
      'Effofecks'; 'Derivate '; 'Electron ';
      'Wise     '; 'Muscle   '; 'Biceps   ' ];
   s_names = ['English    '; 'Biology    '; 'Histo-Geo  '; 
      'Mathematics'; 'Physics    '; 'Philosophy '; 
      'Sport      ' ];
   subject = [ 1  2  3  4  4  5  6  7  7]';
   d_names = ['Mon'; 'Tue'; 'Wed'; 'Thu'; 'Fri' ];
   l_times = ['0800-1000'; '1015-1215'; '1400-1600'; '1615-1815'; ];
   t_nr   = length(t_names);
   s_nr   = 20 ; % number of slots
   temp   = reshape(Result.x_k,t_nr*2,s_nr);
   class1 = temp(1:t_nr,:);
   class2 = temp(t_nr+1:end,:);
   classes= [class1 class2];
   for c = 1:2,
      this_class = class1;
      if c == 2,
         this_class = class2;
      end
      disp(' ')
      disp(['TIMETABLE FOR CLASS ' num2str(c)])
      counter = 0;
      for i = 1:length(d_names),
         day = d_names(i,:);
         disp(['--' day '--'])
         for j = 1:size(l_times,1),
            time = l_times(j,:);
            counter = counter + 1;
            if sum(this_class(:,counter)) == 1, 
               teacher = find(this_class(:,counter));
               disp(['   ' time ' ' s_names(subject(teacher),:) ...
                     ' (' t_names(teacher,:) ')'])
            end
         end
      end
   end 
end

% MODIFICATION LOG
%
% 051205 med   Created.
% 060117 per   Added documentation.
% 060126 per   Moved disp to end

%% Establishing a College Timetable
%
%% Problem description
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
%
%  +------------+-----------------+--------------+--------------+
%  |Teacher     |Subject          |Lsns for cls 1|Lsns for cls 2|
%  +------------+-----------------+--------------+--------------+
%  |Mr Cheese   |English          |      1       |      1       |
%  |Mrs Insulin |Biology          |      3       |      3       |
%  |Mr Map      |History-Geography|      2       |      2       |
%  |Mr Effofecks|Mathematics      |      0       |      4       |
%  |Mrs Derivate|Mathematics      |      4       |      0       |
%  |Mrs Electron|Physics          |      3       |      3       |
%  |Mr Wise     |Philosophy       |      1       |      1       |
%  |Mr Muscle   |Sport            |      1       |      0       |
%  |Mrs Biceps  |Sport            |      0       |      1       |
%  +------------+-----------------+--------------+--------------+
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
%% Variables
%
%  lessons1/2                 Number of lessons per subject and class
%  subject                    Subject indices
%  slots                      Possible slots
%
%% References
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Dec 5, 2005.   Last modified Apr 11, 2009.

%% Problem setup
lessons1        = [ 1  3  2  0  4  3  1  1  0]';
lessons2        = [ 1  3  2  4  0  3  1  0  1]';
subject         = [ 1  2  3  4  4  5  6  7  7]';
slots           = 4*5;

n1    = length(lessons1);   %teachers
n2    = 2;                  %2
n3    = slots;              %slots

t = tomArrayIdx('t',1:n1); % index over all teachers
c = tomArrayIdx('c',1:n2); % index over all classes
l = tomArrayIdx('l',1:n3); % index over all lessons

teach = tom('teach',n1*n2*n3,1,'int');
teach = tomArray(teach,[n1,n2,n3]);

% All variables are binary
bnds1 = {0 <= teach <= 1};

% Sport has to be on thursday afternoon
bnds2 = {1 <= teach(8,1,15), 1 <= teach(9,2,15)};

% No lessons Monday morning
bnds3 = {teach(t,c,1) <= 0};

% Teacher 4 is away on Monday morning
bnds4 = {teach(4,c,2) <= 0};

% Teacher 2 is away on Wednesdays
bnds5 = {teach({'t',2},c,{'l',9:12}) <= 0};

bnds = {bnds1, bnds2, bnds3, bnds4, bnds5};

% Course constraints
con1 = {sum(teach(t,c,l),l) == [lessons1 lessons2]};

% Teacher constraint, one teacher per slot and class
con2 = {sum(teach(t,c,l),t) <= 1};

% Class constraint, one class per slot
con3 = {sum(teach(t,c,l),c) <= 1};

% At most one 2-hour slot of a subject per day
con4 = cell(5,1);
con5 = cell(5,1);
con6 = cell(5,1);
con7 = cell(5,1);

for weekit=1:5
    wl = tomArrayIdx('l',(weekit-1)*4+1:weekit*4); % index over slots
    %Biology
    con4{weekit} = {sum(teach(2,c,wl),wl) <= 1};
    %History
    con5{weekit} = {sum(teach(3,c,wl),wl) <= 1};
    %Math
    wt = tomArrayIdx('t',4:5);
    con6{weekit} = {sum(teach(wt,c,wl),wl) <= 1};
    %Physics
    con7{weekit} = {sum(teach(6,c,wl),wl) <= 1};
end

con = {con1, con2, con3, con4, con5, con6, con7};

% Objective
mornaftslots = tomArrayIdx('l',[1 4 5 8 9 12 13 16 17 20]);
objective = sum(vec(teach(t,c,mornaftslots)));

options = struct;
options.solver = 'cplex';
options.name   = 'College Time Table';
sol = ezsolve(objective,{bnds, con},[],options);

PriLev = 1;
if PriLev > 0
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
    temp   = reshape(sol.teach,t_nr*2,s_nr);
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
            d_name = d_names(i,:);
            disp(['--' d_name '--'])
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
% 051205 med   Created
% 060117 per   Added documentation
% 060126 per   Moved disp to end
% 090412 med   Converted to tomSym
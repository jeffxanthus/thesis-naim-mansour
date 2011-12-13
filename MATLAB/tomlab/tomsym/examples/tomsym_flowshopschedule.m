%% Flow-Shop Scheduling
%
%% Problem description
% A workshop that produces metal pipes on demand for the automobile
% industry has three machines for bending the pipes, soldering the
% fastenings, and assembling the links. The workshop has to produce
% six pieces, for which the durations of the processing steps (in
% minutes) are given in the following table. Every workpiece first
% goes to bending, then to soldering, and finally to assembly of the
% links. Once started, any operations must be carried out without
% interruption, but the workpieces may wait between the machines.
%
% Processing durations in minutes
%
%  +---------+-+-+-+-+-+-+
%  |Workpiece|1|2|3|4|5|6|
%  +---------+-+-+-+-+-+-+
%  |Bending  |3|6|3|5|5|7|
%  |Soldering|5|4|2|4|4|5|
%  |Assembly |5|2|4|6|3|6|
%  +---------+-+-+-+-+-+-+
%
% Every machine only processes one piece at a time. A workpiece may
% not overtake any other by passing onto the following machine. This
% means that if at the beginning a sequence of the workpieces is
% established, they will be processed on every machine in exactly
% this order. Which is the sequence of workpieces that minimizes the
% total time for completing all pieces?
%
%% Variables
%
%  nummachines                The number of machines
%  proctimes                  Time to process a workpiece in a machine
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
nummachines = 3;
proctimes   = [3 6 3 5 5 7;...
    5 4 2 4 4 5;...
    5 2 4 6 3 6];

n1 = size(proctimes,2); % Number of jobs     (j in JOBS)
n2 = size(proctimes,2); % Number of ranks    (k in RANKS)
n3 = nummachines;       % Number of machines (m in MACH)

rnk = tom('rnk',n1,n2,'int');
empty = tom('empty',n3,n2);
wait = tom('wait',n3,n2);

% All slots are integers.
bnds1 = {0 <= rnk <= 1};

% Assignment constraint one job per rank.
con1 = {sum(rnk,1) == 1};

% Assignment constraint one rank per job.
con2 = {sum(rnk,2) == 1};

% Relationship between the end of job rank k on machine m and start of job
% on machines m+1

duration = tom('duration',n3,n2);
durcon = {};
for m=1:n3
    for k=1:n2
        durcon{(m-1)*n2+k} = {duration(m,k) == sum(proctimes(m,:)*rnk(:,k))};
    end
end

con3 = {empty(1:end-1,1:end-1) + duration(1:end-1,2:end) + ...
    wait(1:end-1,2:end) == wait(1:end-1,1:end-1) + ...
    duration(2:end,1:end-1) + empty(2:end,1:end-1)};

% Empty and Wait are zeros when starting, set in bounds
bnds2 = {empty(:,1:end-1) >= 0};
bnds3 = {wait(1:end-1,:) >= 0};
bnds4 = {empty(1,1:end-1) == 0};
bnds5 = {wait(1:end-1,1) == 0};

% Objective
objective = sum(sum(proctimes(1:end-1,:)*rnk(1:end,1))) + sum(sum(empty(end,1:end-1)));

constraints = {bnds1, bnds2, bnds3, bnds4, bnds5, con1, con2, con3, durcon};
options = struct;
options.solver = 'cplex';
options.name   = 'Flow Shop Scheduling';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    p      = size(proctimes,2);         % number of pieces
    m      = nummachines;               % number of machines
    order = round(sol.rnk);
    wait1  = round(sol.wait(1,:));      % wait before machine 1
    wait2  = round(sol.wait(2,:));      % wait before machine 2
    wait3  = round(sol.wait(3,:));      % wait before machine 3
    proctimes   = [3 6 3 5 5 7;...      % the processing times
        5 4 2 4 4 5;...
        5 2 4 6 3 6];
    sequence = [] ;                     % the sequence of pieces
    for i = 1:p,
        sequence = [sequence find(order(:,i))]; % insert piece by piece
    end
    disp(['A best sequence is: ' num2str(sequence)])

    mt1 = [];   % empty machine-times for machine 1
    mt2 = [];   % empty machine-times for machine 2
    mt3 = [];   % empty machine-times for machine 3
    t11 = 0;    % timepoint for entering machine 1
    t12 = 0;    % timepoint for exiting  machine 1
    t21 = 0;    % timepoint for entering machine 2
    t22 = 0;    % timepoint for exiting  machine 2
    t31 = 0;    % timepoint for entering machine 3
    t32 = 0;    % timepoint for exiting  machine 3
    first2 = 1; % the first piece in machine 2
    first3 = 1; % the first piece in machine 3

    for i = 1:length(sequence),
        id  = sequence(i);

        % times in and out of machine 1
        t11  = t12 + wait1(id);
        t12  = t11 + proctimes(1,id);

        % times in and out of machine 2
        if first2 == 1,
            t21 = t12 + wait2(id);
            first2 = 0 ;
        else
            t21 = t22 + wait2(id);
        end
        t22  = t21 + proctimes(2,id);

        % times in and out of machine 3
        if first3 == 1,
            t31 = t22 + wait3(id);
            first3 = 0 ;
        else
            t31 = t32 + wait3(id);
        end
        t32  = t31 + proctimes(3,id);

        % times in and out of machines
        mt1 = [mt1; t11  t12];
        mt2 = [mt2; t21  t22];
        mt3 = [mt3; t31  t32];
    end

    mt = [mt1 mt2 mt3];

    disp('FLOW FOR THE PIECES')
    for j = 1:p,
        disp(['piece ' num2str(sequence(j)) ' has this flow' ])
        for k = 1:m,
            disp(['   machine ' num2str(k) ': ' num2str(mt(j,(k-1)*2+1)) '-' num2str(mt(j,(k-1)*2+2)) ])
        end
    end
    disp('FLOW FOR THE MACHINES')
    for k = 1:m,
        disp(['machine ' num2str(k) ' has this flow' ])
        for l = 1:p,
            disp(['   piece ' num2str(sequence(l)) ': ' ...
                num2str(mt(l,((k-1)*2+1))) '-' num2str(mt(l,((k-1)*2+2)))])
        end
    end
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060111 per   Added documentation.
% 060124 per   Interpretation of results upgraded.
% 060126 per   Moved disp to end
% 090308 med   Converted to tomSym